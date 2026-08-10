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

        sequence = str(record.seq).strip()

        if protein_id in uniprot_ids and sequence:
            rows.append({
                "Protein_ID": protein_id,
                "IEDB_Sequence": sequence,
            })

    sequence_df = pd.DataFrame(
        rows,
        columns=["Protein_ID", "IEDB_Sequence"],
    )

    if sequence_df.empty:
        return sequence_df

    sequence_df["IEDB_Sequence"] = (
        sequence_df["IEDB_Sequence"]
        .astype("string")
        .str.strip()
    )

    sequence_counts = (
        sequence_df
        .dropna(subset=["Protein_ID", "IEDB_Sequence"])
        .groupby("Protein_ID")["IEDB_Sequence"]
        .nunique()
    )
    conflicting_ids = sequence_counts[sequence_counts > 1].index.tolist()

    if conflicting_ids:
        examples = ", ".join(conflicting_ids[:20])
        raise ValueError(
            "Conflicting IEDB sequences were found for the same Protein_ID "
            f"in {fasta_path}. Example IDs: {examples}"
        )

    return sequence_df.drop_duplicates(subset=["Protein_ID"], keep="first")


def ensure_iedb_fasta_exists(fasta_path, protein_ids):
    protein_ids = set(protein_ids)
    fasta_path.parent.mkdir(parents=True, exist_ok=True)

    existing_ids = set()
    if fasta_path.exists():
        for record in SeqIO.parse(fasta_path, "fasta"):
            protein_id = (
                record.id.split("|")[1]
                if "|" in record.id
                else record.id.split()[0]
            )
            protein_id = clean_id(protein_id)
            if protein_id and str(record.seq).strip():
                existing_ids.add(protein_id)

        print(
            f"Existing IEDB FASTA contains {len(existing_ids)} unique "
            "protein IDs."
        )

    missing_ids = sorted(protein_ids - existing_ids)

    if not missing_ids:
        print(f"IEDB FASTA already covers all required proteins: {fasta_path}")
        return []

    print(
        "Fetching "
        f"{len(missing_ids)} missing IEDB source protein sequence(s) "
        "from UniProt..."
    )

    failed = []

    mode = "a" if fasta_path.exists() else "w"
    with open(fasta_path, mode, encoding="utf-8") as out_f:
        if mode == "a" and fasta_path.stat().st_size > 0:
            out_f.write("\n")

        for protein_id in tqdm(missing_ids, desc="Fetching UniProt FASTA"):
            fasta_text = fetch_uniprot_fasta(protein_id)

            if not fasta_text:
                failed.append(protein_id)
                continue

            out_f.write(fasta_text.strip() + "\n")

    print(f"Fetched FASTA sequences for {len(missing_ids) - len(failed)} proteins.")

    if failed:
        print(f"Warning: failed to fetch {len(failed)} proteins from UniProt.")

    return failed


def validate_required_columns(df, required_columns, table_name):
    missing_columns = sorted(set(required_columns) - set(df.columns))

    if missing_columns:
        raise ValueError(
            f"{table_name} is missing required columns: "
            f"{', '.join(missing_columns)}"
        )


def build_epitope_source_map(match_df):
    source_map = match_df[["IEDB_Protein_ID", "Epitope_Source"]].copy()
    source_map["Epitope_Source"] = (
        source_map["Epitope_Source"]
        .astype("string")
        .str.strip()
        .replace("", pd.NA)
    )

    missing_source_ids = (
        source_map.loc[
            source_map["Epitope_Source"].isna(),
            "IEDB_Protein_ID",
        ]
        .dropna()
        .unique()
        .tolist()
    )

    if missing_source_ids:
        examples = ", ".join(sorted(missing_source_ids)[:20])
        raise ValueError(
            "Missing Epitope_Source labels for IEDB protein IDs. "
            f"Example IDs: {examples}"
        )

    source_counts = (
        source_map
        .dropna(subset=["IEDB_Protein_ID", "Epitope_Source"])
        .groupby("IEDB_Protein_ID")["Epitope_Source"]
        .nunique()
    )
    conflicting_ids = source_counts[source_counts > 1].index.tolist()

    if conflicting_ids:
        conflicts = (
            source_map[source_map["IEDB_Protein_ID"].isin(conflicting_ids)]
            .drop_duplicates()
            .sort_values(["IEDB_Protein_ID", "Epitope_Source"])
        )
        print("Conflicting Epitope_Source labels:")
        print(conflicts.head(50).to_string(index=False))
        raise ValueError(
            "Each IEDB_Protein_ID must map to exactly one "
            "Epitope_Source label."
        )

    return source_map.drop_duplicates(subset=["IEDB_Protein_ID"])


def build_pathogen_sequence_table(pathogen_data):
    pathogen_data = pathogen_data.copy()
    pathogen_data["Protein_ID"] = pathogen_data["Protein_ID"].map(clean_id)

    missing_id_rows = pathogen_data["Protein_ID"].isna().sum()
    if missing_id_rows:
        raise ValueError(
            f"protein_sequences.csv contains {missing_id_rows} row(s) "
            "with missing Protein_ID."
        )

    duplicated_ids = pathogen_data.loc[
        pathogen_data["Protein_ID"].duplicated(keep=False),
        "Protein_ID",
    ].unique()

    if len(duplicated_ids) > 0:
        examples = ", ".join(sorted(duplicated_ids)[:20])
        raise ValueError(
            "protein_sequences.csv must contain exactly one row per "
            f"Protein_ID. Duplicate IDs include: {examples}"
        )

    pathogen_data["Sequence"] = (
        pathogen_data["Sequence"]
        .astype("string")
        .str.strip()
        .replace("", pd.NA)
    )

    missing_sequence_ids = pathogen_data.loc[
        pathogen_data["Sequence"].isna(),
        "Protein_ID",
    ].tolist()

    if missing_sequence_ids:
        examples = ", ".join(sorted(missing_sequence_ids)[:20])
        raise ValueError(
            "protein_sequences.csv contains proteins without a sequence. "
            f"Example IDs: {examples}"
        )

    return (
        pathogen_data[
            ["Protein_ID", "Sequence", "Annotation", "Gene_name"]
        ]
        .rename(columns={
            "Protein_ID": "Pathogen_Protein_ID",
            "Sequence": "Pathogen_Sequence",
            "Annotation": "Pathogen_Annotation",
            "Gene_name": "Pathogen_Gene_Name",
        })
    )


def empty_alignment_metrics(metric_cols):
    return pd.Series({col: None for col in metric_cols})


def compute_alignment_metrics(seq1, seq2, aligner, metric_cols):
    if pd.isna(seq1) or pd.isna(seq2):
        return empty_alignment_metrics(metric_cols)

    seq1 = str(seq1).strip()
    seq2 = str(seq2).strip()

    if not seq1 or not seq2:
        return empty_alignment_metrics(metric_cols)

    try:
        best = aligner.align(seq1, seq2)[0]
    except Exception:
        return empty_alignment_metrics(metric_cols)

    # Let Biopython count the alignment components
    counts = best.counts()

    matches = counts.identities
    mismatches = counts.mismatches
    gap_count = counts.gaps
    alignment_length = best.length

    if alignment_length == (matches + mismatches + gap_count):
        pass  # Valid alignment
    else:
        print("Warning: alignment length mismatch.")
        print(f"Seq1: {seq1}")
        print(f"Seq2: {seq2}")
        print(f"Alignment:\n{best}")
        print(f"Counts: {counts}")
        return empty_alignment_metrics(metric_cols)

    # Gap-inclusive percentage identity
    percent_identity = (
        matches / (alignment_length) * 100
        if alignment_length > 0
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
        "Normalized_Alignment_Score": normalized_score
    })


# ---------------------- Load data ---------------------- #

print("Loading data...")

match_df = load_csv(DATA_DIR / "proccesed/iedb_match_regions_long.csv")
pathogen_data = load_csv(DATA_DIR / "intermediate/protein_sequences.csv")

validate_required_columns(
    match_df,
    ["IEDB_Protein_ID", "Pathogen_Protein_ID", "Epitope_Source"],
    "iedb_match_regions_long.csv",
)
validate_required_columns(
    pathogen_data,
    ["Protein_ID", "Sequence", "Annotation", "Gene_name"],
    "protein_sequences.csv",
)

print(f"Long match rows: {len(match_df)}")


# ---------------------- Normalize protein IDs ---------------------- #

match_df["Pathogen_Protein_ID"] = match_df["Pathogen_Protein_ID"].map(clean_id)
match_df["IEDB_Protein_ID"] = match_df["IEDB_Protein_ID"].map(clean_id)

match_df = match_df.dropna(
    subset=["Pathogen_Protein_ID", "IEDB_Protein_ID"]
).copy()

print(f"Rows after dropping missing protein IDs: {len(match_df)}")

protein_pairs = (
    match_df[["IEDB_Protein_ID", "Pathogen_Protein_ID"]]
    .drop_duplicates()
    .sort_values(["IEDB_Protein_ID", "Pathogen_Protein_ID"])
    .reset_index(drop=True)
)

epitope_source_map = build_epitope_source_map(match_df)
pathogen_sequences = build_pathogen_sequence_table(pathogen_data)

print(f"Unique protein-protein pairs for alignment: {len(protein_pairs)}")


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


# ---------------------- Prepare sequence table ---------------------- #

iedb_seq = uniprot_seq_df.rename(columns={
    "Protein_ID": "IEDB_Protein_ID",
})

print(f"Unique IEDB proteins with sequences: {len(iedb_seq)}")


# ---------------------- Build protein-pair alignment table ---------------------- #

full_align_analysis = (
    protein_pairs
    .merge(
        epitope_source_map,
        how="left",
        on="IEDB_Protein_ID",
        validate="many_to_one",
    )
    .merge(
        pathogen_sequences,
        how="left",
        on="Pathogen_Protein_ID",
        validate="many_to_one",
    )
    .merge(
        iedb_seq,
        how="left",
        on="IEDB_Protein_ID",
        validate="many_to_one",
    )
)

if len(full_align_analysis) != len(protein_pairs):
    raise RuntimeError(
        "Protein-pair row count changed while attaching sequence metadata."
    )

missing_pathogen_ids = (
    full_align_analysis.loc[
        full_align_analysis["Pathogen_Sequence"].isna(),
        "Pathogen_Protein_ID",
    ]
    .drop_duplicates()
    .tolist()
)

if missing_pathogen_ids:
    examples = ", ".join(sorted(missing_pathogen_ids)[:20])
    raise ValueError(
        "Matched pathogen proteins were not found in protein_sequences.csv. "
        f"Example IDs: {examples}"
    )

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
    "Normalized_Alignment_Score"
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
        Mean_Alignment_Length=("Alignment_Length", "mean"),
        Mean_Gap_Count=("Gap_Count", "mean"),
        N_pathogen_proteins=("Pathogen_Protein_ID", "nunique"),
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

sc = plt.scatter(
    summary_df["Mean_Identity"],
    summary_df["Mean_Normalized_Alignment_Score"],
    c=summary_df["N_pathogen_proteins"],
    s=60,
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
    s=80,
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

plt.tight_layout()

plot_path = DATA_DIR / "proccesed/global_alignment_summary_plot.png"
plt.savefig(plot_path, dpi=300, bbox_inches="tight")
print(f"Saved plot to: {plot_path}")

plt.show()