import pandas as pd
import matplotlib.pyplot as plt

from adjustText import adjust_text
from Bio import Align, SeqIO
from Bio.Align import substitution_matrices
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

    for record in tqdm(
        SeqIO.parse(fasta_path, "fasta"),
        desc="Loading IEDB sequences"
    ):
        protein_id = (
            record.id.split("|")[1]
            if "|" in record.id
            else record.id.split()[0]
        )

        protein_id = clean_id(protein_id)

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

    with open(fasta_path, "w", encoding="utf-8") as out_f:
        for protein_id in tqdm(protein_ids, desc="Fetching UniProt FASTA"):
            fasta_text = fetch_uniprot_fasta(protein_id)

            if fasta_text is None:
                failed.append(protein_id)
                continue

            out_f.write(fasta_text.strip() + "\n")

    print(f"Fetched FASTA sequences for {len(protein_ids) - len(failed)} proteins.")

    if failed:
        print(f"Warning: failed to fetch {len(failed)} proteins from UniProt.")


def empty_alignment_metrics(metric_cols):
    return pd.Series({col: None for col in metric_cols})


def compute_alignment_metrics(seq1, seq2, aligner, metric_cols):
    if pd.isna(seq1) or pd.isna(seq2):
        return empty_alignment_metrics(metric_cols)

    seq1, seq2 = str(seq1).strip(), str(seq2).strip()

    if not seq1 or not seq2:
        return empty_alignment_metrics(metric_cols)

    try:
        best = aligner.align(seq1, seq2)[0]
    except Exception:
        return empty_alignment_metrics(metric_cols)

    matches = 0
    mismatches = 0
    gap_count = 0
    alignment_length = 0

    query_positions = set()
    subject_positions = set()

    coords = best.coordinates

    for i in range(coords.shape[1] - 1):
        start1, end1 = coords[0, i], coords[0, i + 1]
        start2, end2 = coords[1, i], coords[1, i + 1]

        len1 = end1 - start1
        len2 = end2 - start2

        if len1 > 0 and len2 > 0:
            segment1 = seq1[start1:end1]
            segment2 = seq2[start2:end2]

            query_positions.update(range(start1, end1))
            subject_positions.update(range(start2, end2))

            alignment_length += len(segment1)

            for a, b in zip(segment1, segment2):
                if a == b:
                    matches += 1
                else:
                    mismatches += 1

        else:
            gap_len = max(len1, len2)
            gap_count += gap_len
            alignment_length += gap_len

    percent_identity = (
        (matches / alignment_length) * 100
        if alignment_length
        else 0.0
    )

    query_coverage = (
        (len(query_positions) / len(seq1)) * 100
        if len(seq1)
        else 0.0
    )

    subject_coverage = (
        (len(subject_positions) / len(seq2)) * 100
        if len(seq2)
        else 0.0
    )

    raw_score = best.score
    smallest_sequence_length = min(len(seq1), len(seq2))

    normalized_score = (
        raw_score / smallest_sequence_length
        if smallest_sequence_length > 0
        else None
    )

    return pd.Series({
        "Percent_Identity": percent_identity,
        "Alignment_Length": alignment_length,
        "Matches": matches,
        "Mismatches": mismatches,
        "Gap_Count": gap_count,
        "Alignment_Score": raw_score,
        "Normalized_Alignment_Score": normalized_score,
        "Query_Coverage": query_coverage,
        "Subject_Coverage": subject_coverage,
    })


# ---------------------- Load data ---------------------- #

print("Loading data...")

match_df = load_csv(DATA_DIR / "proccesed/iedb_match_regions_long.csv")
pathogen_data = load_csv(DATA_DIR / "intermediate/protein_metadata.csv")

print(f"Long match rows: {len(match_df)}")


# ---------------------- Normalize protein IDs ---------------------- #

match_df["Pathogen_Protein_ID"] = match_df["Pathogen_Protein_ID"].map(clean_id)
match_df["IEDB_Protein_ID"] = match_df["IEDB_Protein_ID"].map(clean_id)
pathogen_data["Protein_ID"] = pathogen_data["Protein_ID"].map(clean_id)

match_df = match_df.dropna(
    subset=["Pathogen_Protein_ID", "IEDB_Protein_ID"]
).copy()

pathogen_data = pathogen_data.dropna(
    subset=["Protein_ID"]
).copy()

print(f"Rows after dropping missing protein IDs: {len(match_df)}")


# ---------------------- Load/fetch IEDB source protein FASTA ---------------------- #

fasta_path = DATA_DIR / "intermediate" / "matched_iedb_source_proteins.fasta"

matched_iedb_ids = (
    match_df["IEDB_Protein_ID"]
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

    save_csv(
        pd.DataFrame({"IEDB_Protein_ID": sorted(missing_iedb_ids)}),
        DATA_DIR / "proccesed/missing_iedb_ids_for_alignment.csv"
    )


# ---------------------- Prepare merge tables ---------------------- #

pathogen_meta = (
    pathogen_data[
        [
            "Protein_ID",
            "Genus_species",
            "Sequence",
            "Strain",
        ]
    ]
    .drop_duplicates(subset=["Protein_ID"])
    .rename(columns={
        "Protein_ID": "Pathogen_Protein_ID_merge",
        "Genus_species": "Organism_Source",
        "Sequence": "Pathogen_Sequence",
    })
)

iedb_seq = (
    uniprot_seq_df
    .drop_duplicates(subset=["Protein_ID"])
    .rename(columns={
        "Protein_ID": "IEDB_Protein_ID_merge",
    })
)

print(f"Unique IEDB proteins with sequences: {len(iedb_seq)}")


# ---------------------- Build full alignment table ---------------------- #

traceback_cols = [
    "IEDB_Protein_ID",
    "Pathogen_Protein_ID",
    "Epitope_Source",
    "Disease",
    "Disease_stage",
    "Response_measured",
    "Effector_cell",
    "Pathogen_Organism",
    "Pathogen_Scientific_name",
    "Pathogen_Annotation",
    "Pathogen_Gene_Name",
    "IEDB_Region_ID",
    "IEDB_Region_Start",
    "IEDB_Region_End",
    "IEDB_Region_Length",
]

traceback_cols = [
    col for col in traceback_cols
    if col in match_df.columns
]

protein_pair_traceback = (
    match_df[traceback_cols]
    .drop_duplicates()
)

full_align_analysis = (
    protein_pair_traceback
    .drop_duplicates(subset=["IEDB_Protein_ID", "Pathogen_Protein_ID"])
    .merge(
        pathogen_meta,
        how="left",
        left_on="Pathogen_Protein_ID",
        right_on="Pathogen_Protein_ID_merge",
    )
    .merge(
        iedb_seq,
        how="left",
        left_on="IEDB_Protein_ID",
        right_on="IEDB_Protein_ID_merge",
    )
    .drop(columns=["Pathogen_Protein_ID_merge", "IEDB_Protein_ID_merge"])
)

print(f"Unique protein-protein pairs for alignment: {len(full_align_analysis)}")

unique_pathogen_proteins = full_align_analysis["Pathogen_Protein_ID"].nunique()
print(f"Unique pathogen proteins: {unique_pathogen_proteins}")


# ---------------------- Alignment setup ---------------------- #

aligner = Align.PairwiseAligner()
aligner.mode = "global"
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

aligner.open_gap_score = -5
aligner.extend_gap_score = -1

print("Open:", aligner.open_gap_score)
print("Extend:", aligner.extend_gap_score)

metric_cols = [
    "Percent_Identity",
    "Alignment_Length",
    "Matches",
    "Mismatches",
    "Gap_Count",
    "Alignment_Score",
    "Normalized_Alignment_Score",
    "Query_Coverage",
    "Subject_Coverage",
]


# ---------------------- Compute global alignment metrics ---------------------- #

print("Computing global alignment metrics...")

metrics = full_align_analysis.progress_apply(
    lambda row: compute_alignment_metrics(
        row["IEDB_Sequence"],
        row["Pathogen_Sequence"],
        aligner,
        metric_cols,
    ),
    axis=1,
)

full_align_analysis = pd.concat([full_align_analysis, metrics], axis=1)

alignment_output_path = DATA_DIR / "proccesed/global_alignment_analysis.csv"
save_csv(full_align_analysis, alignment_output_path)

print(f"Saved full alignment analysis to: {alignment_output_path}")


# ---------------------- Summary table ---------------------- #

summary_df = (
    full_align_analysis
    .dropna(subset=["Percent_Identity", "Normalized_Alignment_Score"])
    .groupby(["IEDB_Protein_ID", "Epitope_Source"])
    .agg(
        Mean_Identity=("Percent_Identity", "mean"),
        SD_Identity=("Percent_Identity", "std"),
        Mean_Alignment_Score=("Alignment_Score", "mean"),
        SD_Alignment_Score=("Alignment_Score", "std"),
        Mean_Normalized_Alignment_Score=("Normalized_Alignment_Score", "mean"),
        SD_Normalized_Alignment_Score=("Normalized_Alignment_Score", "std"),
        Mean_Query_Coverage=("Query_Coverage", "mean"),
        Mean_Subject_Coverage=("Subject_Coverage", "mean"),
        Mean_Alignment_Length=("Alignment_Length", "mean"),
        Mean_Gap_Count=("Gap_Count", "mean"),
        N_pathogen_proteins=("Pathogen_Protein_ID", "nunique"),
        N_protein_pairs=("IEDB_Protein_ID", "count"),
    )
    .reset_index()
    .fillna(0)
    .sort_values("N_pathogen_proteins", ascending=False)
    .round(2)
)

summary_output_path = DATA_DIR / "proccesed/global_alignment_summary.csv"
save_csv(summary_df, summary_output_path)

print(f"Saved alignment summary to: {summary_output_path}")


# ---------------------- Scatter plot: global protein similarity ---------------------- #

plt.figure(figsize=(9, 7))

size_scale = 3

sc = plt.scatter(
    summary_df["Mean_Identity"],
    summary_df["Mean_Normalized_Alignment_Score"],
    c=summary_df["N_pathogen_proteins"],
    s=summary_df["Mean_Query_Coverage"] * size_scale,
    norm=LogNorm(),
    alpha=0.75,
)

plt.errorbar(
    summary_df["Mean_Identity"],
    summary_df["Mean_Normalized_Alignment_Score"],
    xerr=summary_df["SD_Identity"],
    yerr=summary_df["SD_Normalized_Alignment_Score"],
    fmt="none",
    alpha=0.25,
    capsize=2,
)

highlight_ids = [
    # Already selected / major proteins
    "P10809",  # 60 kDa heat shock protein, mitochondrial
    "P0DMV8",  # Heat shock 70 kDa protein 1A
    "P11021",  # Endoplasmic reticulum chaperone BiP
    "P01308",  # Insulin
    "P62805",  # Histone H4
    "P38646",  # Stress-70 protein, mitochondrial
    "Q71DI3",  # Histone H3.2
    "P62807",  # Histone H2B type 1-C/E/F/G/I
    "P11142",  # Heat shock cognate 71 kDa protein
    "P55087",  # Aquaporin-4
    "Q05329",  # Glutamate decarboxylase 2
    "P06733",  # Alpha-enolase

    # Additional selected proteins
    "P02671",  # Fibrinogen alpha chain
    "Q8IWU4",  # Proton-coupled zinc antiporter SLC30A8
    "P02458",  # Collagen alpha-1(II) chain
    "Q01955",  # Collagen alpha-3(IV) chain
    "P07202",  # Thyroid peroxidase
]

label_df = summary_df[summary_df["IEDB_Protein_ID"].isin(highlight_ids)]

plt.scatter(
    label_df["Mean_Identity"],
    label_df["Mean_Normalized_Alignment_Score"],
    s=label_df["Mean_Query_Coverage"] * size_scale * 1.3,
    facecolors="none",
    edgecolors="black",
    linewidths=1.2,
    alpha=0.9,
)

texts = [
    plt.text(
        row["Mean_Identity"],
        row["Mean_Normalized_Alignment_Score"],
        row["Epitope_Source"],
        fontsize=8,
        bbox=dict(
            boxstyle="round,pad=0.2",
            facecolor="white",
            edgecolor="none",
            alpha=0.7,
        ),
    )
    for _, row in label_df.iterrows()
]

adjust_text(
    texts,
    x=label_df["Mean_Identity"].to_numpy(),
    y=label_df["Mean_Normalized_Alignment_Score"].to_numpy(),
    arrowprops=dict(
        arrowstyle="-",
        color="gray",
        lw=0.6,
        alpha=0.45,
        shrinkA=3,
        shrinkB=3,
    ),
    force_text=(0.8, 0.8),
    force_points=(0.4, 0.4),
    expand_points=(1.3, 1.3),
    expand_text=(1.2, 1.2),
)

plt.xlabel("Mean global percent identity (%)")
plt.ylabel("Mean normalized global alignment score")

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
    frameon=True,
)

plt.tight_layout()

plot_path = DATA_DIR / "proccesed/global_alignment_summary_plot.png"
plt.savefig(plot_path, dpi=300, bbox_inches="tight")
print(f"Saved plot to: {plot_path}")

plt.show()