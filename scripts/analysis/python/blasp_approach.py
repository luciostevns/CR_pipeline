import shutil
import subprocess

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from adjustText import adjust_text
from Bio import SeqIO
from tqdm import tqdm
from matplotlib.colors import LogNorm

from crossreactivity.io import load_csv, save_csv, DATA_DIR, clean_id
from crossreactivity.uniprot import fetch_uniprot_fasta


tqdm.pandas()


# ---------------------- Helper functions ---------------------- #

def load_iedb_sequences_from_fasta(fasta_path, uniprot_ids):
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    uniprot_ids = set(uniprot_ids)
    rows = []

    for record in tqdm(SeqIO.parse(fasta_path, "fasta"), desc="Loading IEDB sequences"):
        protein_id = record.id.split("|")[1] if "|" in record.id else record.id.split()[0]

        if protein_id in uniprot_ids:
            rows.append({
                "Protein_ID": protein_id,
                "IEDB_Sequence": str(record.seq)
            })

    return pd.DataFrame(rows)


def ensure_iedb_fasta_exists(fasta_path, protein_ids):
    if fasta_path.exists():
        print(f"Using existing FASTA: {fasta_path}")
        return

    print("Fetching matched IEDB source protein sequences from UniProt...")
    fasta_path.parent.mkdir(parents=True, exist_ok=True)

    failed = []

    with open(fasta_path, "w") as out_f:
        for protein_id in tqdm(protein_ids, desc="Fetching UniProt FASTA"):
            fasta_text = fetch_uniprot_fasta(protein_id)

            if fasta_text is None:
                failed.append(protein_id)
                continue

            out_f.write(fasta_text.strip() + "\n")

    print(f"Fetched FASTA sequences for {len(protein_ids) - len(failed)} proteins.")

    if failed:
        print(f"Warning: failed to fetch {len(failed)} proteins from UniProt.")


def check_blast_installed():
    missing = [cmd for cmd in ["makeblastdb", "blastp"] if shutil.which(cmd) is None]

    if missing:
        raise RuntimeError(
            f"Missing BLAST+ command(s): {', '.join(missing)}. "
            "Install BLAST+ and make sure it is available in PATH."
        )


def write_fasta(df, id_col, seq_col, output_path):
    output_path.parent.mkdir(parents=True, exist_ok=True)

    clean_df = (
        df[[id_col, seq_col]]
        .dropna()
        .drop_duplicates(subset=[id_col])
    )

    with open(output_path, "w") as f:
        for _, row in clean_df.iterrows():
            protein_id = str(row[id_col]).strip()
            sequence = str(row[seq_col]).replace(" ", "").replace("\n", "").strip()

            if not protein_id or not sequence:
                continue

            f.write(f">{protein_id}\n")

            for i in range(0, len(sequence), 80):
                f.write(sequence[i:i + 80] + "\n")

    print(f"Wrote {len(clean_df)} sequences to {output_path}")


def run_blastp(query_fasta, subject_fasta, db_prefix, output_tsv, num_threads=4):
    check_blast_installed()

    db_prefix.parent.mkdir(parents=True, exist_ok=True)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    print("Building BLAST database...")

    subprocess.run(
        [
            "makeblastdb",
            "-in", str(subject_fasta),
            "-dbtype", "prot",
            "-out", str(db_prefix),
            "-parse_seqids"
        ],
        check=True
    )

    outfmt = (
        "6 qseqid sseqid pident length mismatch gapopen "
        "qstart qend sstart send evalue bitscore qlen slen qcovs"
    )

    print("Running BLASTp...")

    subprocess.run(
        [
            "blastp",
            "-query", str(query_fasta),
            "-db", str(db_prefix),
            "-out", str(output_tsv),
            "-outfmt", outfmt,
            "-evalue", "10",
            "-max_hsps", "1",
            "-num_threads", str(num_threads)
        ],
        check=True
    )


def load_blast_results(blast_tsv):
    cols = [
        "IEDB_Protein_ID",
        "Pathogen_Protein_ID",
        "Percent_Identity",
        "Alignment_Length",
        "Mismatches",
        "Gap_Openings",
        "Query_Start",
        "Query_End",
        "Subject_Start",
        "Subject_End",
        "Evalue",
        "Bit_Score",
        "Query_Length",
        "Subject_Length",
        "Query_Coverage"
    ]

    if not blast_tsv.exists() or blast_tsv.stat().st_size == 0:
        print("Warning: BLASTp output is empty.")
        return pd.DataFrame(columns=cols + ["Subject_Coverage"])

    blast_df = pd.read_csv(blast_tsv, sep="\t", names=cols)

    blast_df["Subject_Coverage"] = (
        (blast_df["Subject_End"] - blast_df["Subject_Start"]).abs() + 1
    ) / blast_df["Subject_Length"] * 100

    return blast_df


# ---------------------- Load data ---------------------- #

print("Loading data...")

perfect_match = load_csv(DATA_DIR / "proccesed/perfect_matches_2_0.csv")
pathogen_data = load_csv(DATA_DIR / "intermediate/protein_metadata.csv")


# ---------------------- Filter matched rows ---------------------- #

perfect_match["Matched"] = perfect_match["Matched"].astype(str).str.lower() == "true"
perfect_match = perfect_match[perfect_match["Matched"]].copy()

print(f"Matched rows: {len(perfect_match)}")


# ---------------------- Normalize protein IDs ---------------------- #

perfect_match["Pathogen_Protein_ID"] = perfect_match["Pathogen_Protein_ID"].map(clean_id)
perfect_match["IEDB_Protein_ID"] = perfect_match["IEDB_Protein_ID"].map(clean_id)
pathogen_data["Protein_ID"] = pathogen_data["Protein_ID"].map(clean_id)

perfect_match = perfect_match.dropna(
    subset=["Pathogen_Protein_ID", "IEDB_Protein_ID"]
).copy()

pathogen_data = pathogen_data.dropna(
    subset=["Protein_ID"]
).copy()


# ---------------------- Load/fetch IEDB source protein FASTA ---------------------- #

fasta_path = DATA_DIR / "intermediate" / "matched_iedb_source_proteins.fasta"

matched_iedb_ids = (
    perfect_match["IEDB_Protein_ID"]
    .dropna()
    .unique()
    .tolist()
)

ensure_iedb_fasta_exists(
    fasta_path=fasta_path,
    protein_ids=matched_iedb_ids
)

print("Loading IEDB source protein sequences from FASTA...")

uniprot_seq_df = load_iedb_sequences_from_fasta(
    fasta_path=fasta_path,
    uniprot_ids=matched_iedb_ids
)

print(f"Loaded sequences for {len(uniprot_seq_df)} IEDB proteins.")

missing_iedb_ids = set(matched_iedb_ids) - set(uniprot_seq_df["Protein_ID"])

if missing_iedb_ids:
    print(
        f"Warning: {len(missing_iedb_ids)} matched IEDB proteins "
        "were not found in the FASTA."
    )


# ---------------------- Prepare merge tables ---------------------- #

pathogen_meta = (
    pathogen_data[["Protein_ID", "Genus_species", "Sequence", "Strain"]]
    .drop_duplicates(subset=["Protein_ID"])
    .rename(columns={
        "Protein_ID": "Pathogen_Protein_ID_merge",
        "Genus_species": "Organism_Source",
        "Sequence": "Pathogen_Sequence"
    })
)

iedb_seq = (
    uniprot_seq_df
    .drop_duplicates(subset=["Protein_ID"])
    .rename(columns={
        "Protein_ID": "IEDB_Protein_ID_merge"
    })
)

print(f"Unique IEDB proteins with sequences: {len(iedb_seq)}")


# ---------------------- Build full protein-pair table ---------------------- #

full_align_analysis = (
    perfect_match
    .drop_duplicates(subset=["IEDB_Protein_ID", "Pathogen_Protein_ID"])
    .merge(
        pathogen_meta,
        how="left",
        left_on="Pathogen_Protein_ID",
        right_on="Pathogen_Protein_ID_merge"
    )
    .merge(
        iedb_seq,
        how="left",
        left_on="IEDB_Protein_ID",
        right_on="IEDB_Protein_ID_merge"
    )
    .drop(columns=["Pathogen_Protein_ID_merge", "IEDB_Protein_ID_merge"])
)

print(f"Unique IEDB-pathogen protein pairs: {len(full_align_analysis)}")
print(f"Unique pathogen proteins: {full_align_analysis['Pathogen_Protein_ID'].nunique()}")


# ---------------------- Write FASTA files for BLASTp ---------------------- #

blast_dir = DATA_DIR / "intermediate" / "blastp_similarity"
iedb_blast_fasta = blast_dir / "iedb_source_proteins.fasta"
pathogen_blast_fasta = blast_dir / "matched_pathogen_proteins.fasta"
pathogen_db_prefix = blast_dir / "matched_pathogen_db"
blast_output = blast_dir / "iedb_vs_pathogen_blastp.tsv"

iedb_fasta_df = (
    full_align_analysis[["IEDB_Protein_ID", "IEDB_Sequence"]]
    .dropna()
    .drop_duplicates(subset=["IEDB_Protein_ID"])
)

pathogen_fasta_df = (
    full_align_analysis[["Pathogen_Protein_ID", "Pathogen_Sequence"]]
    .dropna()
    .drop_duplicates(subset=["Pathogen_Protein_ID"])
)

write_fasta(
    iedb_fasta_df,
    id_col="IEDB_Protein_ID",
    seq_col="IEDB_Sequence",
    output_path=iedb_blast_fasta
)

write_fasta(
    pathogen_fasta_df,
    id_col="Pathogen_Protein_ID",
    seq_col="Pathogen_Sequence",
    output_path=pathogen_blast_fasta
)


# ---------------------- Run BLASTp ---------------------- #

run_blastp(
    query_fasta=iedb_blast_fasta,
    subject_fasta=pathogen_blast_fasta,
    db_prefix=pathogen_db_prefix,
    output_tsv=blast_output,
    num_threads=4
)

blast_df = load_blast_results(blast_output)

print(f"BLASTp alignments returned: {len(blast_df)}")


# ---------------------- Keep only original matched protein pairs ---------------------- #

pair_table = (
    full_align_analysis[["IEDB_Protein_ID", "Pathogen_Protein_ID", "Epitope_Source"]]
    .drop_duplicates()
)

blast_matched_pairs = pair_table.merge(
    blast_df,
    how="left",
    on=["IEDB_Protein_ID", "Pathogen_Protein_ID"]
)

blast_matched_pairs["Has_BLAST_Hit"] = blast_matched_pairs["Bit_Score"].notna()

print(
    f"Matched protein pairs with BLASTp hit: "
    f"{blast_matched_pairs['Has_BLAST_Hit'].sum()} / {len(blast_matched_pairs)}"
)

save_csv(
    blast_matched_pairs,
    DATA_DIR / "proccesed/iedb_pathogen_blastp_pair_results.csv"
)


# ---------------------- Summary table ---------------------- #

summary_df = (
    blast_matched_pairs
    .groupby(["IEDB_Protein_ID", "Epitope_Source"])
    .agg(
        Mean_Identity=("Percent_Identity", "mean"),
        SD_Identity=("Percent_Identity", "std"),
        Mean_Query_Coverage=("Query_Coverage", "mean"),
        SD_Query_Coverage=("Query_Coverage", "std"),
        Mean_Subject_Coverage=("Subject_Coverage", "mean"),
        Mean_Alignment_Length=("Alignment_Length", "mean"),
        Best_Evalue=("Evalue", "min"),
        Best_Bit_Score=("Bit_Score", "max"),
        Mean_Bit_Score=("Bit_Score", "mean"),
        N_matches=("Pathogen_Protein_ID", "count"),
        N_BLAST_hits=("Has_BLAST_Hit", "sum")
    )
    .reset_index()
)

summary_df["Minus_Log10_Best_Evalue"] = -np.log10(
    summary_df["Best_Evalue"].replace(0, np.nextafter(0, 1))
)

summary_df = (
    summary_df
    .fillna(0)
    .sort_values("N_matches", ascending=False)
    .round(2)
)

save_csv(
    summary_df,
    DATA_DIR / "proccesed/iedb_blastp_similarity_summary.csv"
)

print("Saved BLASTp summary table.")


# ---------------------- Scatter plot: BLASTp similarity ---------------------- #

plot_df = summary_df[summary_df["N_BLAST_hits"] > 0].copy()

plt.figure(figsize=(9, 7))

size_scale = 3

sc = plt.scatter(
    plot_df["Mean_Identity"],
    plot_df["Minus_Log10_Best_Evalue"],
    c=plot_df["N_matches"],
    s=plot_df["Mean_Query_Coverage"] * size_scale,
    norm=LogNorm(),
    alpha=0.75
)

plt.errorbar(
    plot_df["Mean_Identity"],
    plot_df["Minus_Log10_Best_Evalue"],
    xerr=plot_df["SD_Identity"],
    fmt="none",
    alpha=0.25,
    capsize=2
)

highlight_ids = [
    "P10809", "P0DMV8", "P11021", "P01308", "P62805", "P38646",
    "Q71DI3", "P62807", "P11142", "P55087", "Q05329", "P06733"
]

label_df = plot_df[plot_df["IEDB_Protein_ID"].isin(highlight_ids)]

plt.scatter(
    label_df["Mean_Identity"],
    label_df["Minus_Log10_Best_Evalue"],
    s=label_df["Mean_Query_Coverage"] * size_scale * 1.3,
    facecolors="none",
    edgecolors="black",
    linewidths=1.2,
    alpha=0.9
)

texts = [
    plt.text(
        row["Mean_Identity"],
        row["Minus_Log10_Best_Evalue"],
        row["Epitope_Source"],
        fontsize=8,
        bbox=dict(
            boxstyle="round,pad=0.2",
            facecolor="white",
            edgecolor="none",
            alpha=0.7
        )
    )
    for _, row in label_df.iterrows()
]

adjust_text(
    texts,
    x=label_df["Mean_Identity"].to_numpy(),
    y=label_df["Minus_Log10_Best_Evalue"].to_numpy(),
    arrowprops=dict(
        arrowstyle="-|>",
        color="black",
        lw=1.0,
        alpha=0.8,
        shrinkA=4,
        shrinkB=4,
        mutation_scale=10
    ),
    force_text=(0.8, 0.8),
    force_points=(0.4, 0.4),
    expand_points=(1.3, 1.3),
    expand_text=(1.2, 1.2)
)

plt.xlabel("Mean BLASTp percent identity (%)")
plt.ylabel("-log10(best E-value)")

cbar = plt.colorbar(sc)
cbar.set_label("Number of matched pathogen proteins (log scale)")

size_values = [25, 50, 75, 100]
size_handles = [
    plt.scatter([], [], s=value * size_scale, alpha=0.75, label=f"{value}%")
    for value in size_values
]

plt.legend(
    handles=size_handles,
    title="Mean query coverage",
    loc="best",
    frameon=True
)

plt.tight_layout()
plt.show()