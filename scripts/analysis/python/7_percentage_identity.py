import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from adjustText import adjust_text
from Bio import Align
from tqdm import tqdm
from matplotlib.colors import LogNorm
from crossreactivity.io import load_csv, save_csv, DATA_DIR, clean_id

tqdm.pandas()

#Load data
print("Loading data...")

perfect_match = load_csv(DATA_DIR / "proccesed/perfect_matches_2_0.csv")
pathogen_data = load_csv(DATA_DIR / "intermediate/protein_metadata.csv")

# Filter to matched rows only
perfect_match["Matched"] = perfect_match["Matched"] == "True"
perfect_match = perfect_match[perfect_match["Matched"]].copy()

print(f"Matched rows: {len(perfect_match)}")

# Normalize IDs
perfect_match["Pathogen_Protein_ID"] = perfect_match["Pathogen_Protein_ID"].map(clean_id)
perfect_match["IEDB_Protein_ID"] = perfect_match["IEDB_Protein_ID"].map(clean_id)
pathogen_data["Protein_ID"] = pathogen_data["Protein_ID"].map(clean_id)

# Drop rows with missing IDs after cleanup
perfect_match = perfect_match.dropna(subset=["Pathogen_Protein_ID", "IEDB_Protein_ID"]).copy()
pathogen_data = pathogen_data.dropna(subset=["Protein_ID"]).copy()


# Fetch matched IEDB source protein sequences
def fetch_iedb_sequences(uniprot_ids):
    fasta_path = DATA_DIR / "intermediate" / "matched_iedb_source_proteins.fasta"
    
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
    
    rows = []
    
    # Read the FASTA file
    with open(fasta_path, 'r') as f:
        fasta_content = f.read()
    
    # Split by header lines (lines starting with ">")
    entries = []
    current_header = None
    current_seq = []
    
    for line in fasta_content.splitlines():
        if line.startswith(">"):
            if current_header is not None:
                entries.append((current_header, "".join(current_seq)))
            current_header = line[1:].strip()
            current_seq = []
        else:
            current_seq.append(line.strip())
    
    # Don't forget the last entry
    if current_header is not None:
        entries.append((current_header, "".join(current_seq)))
    
    # Match entries to requested uniprot_ids
    for protein_id in tqdm(uniprot_ids, desc="Loading IEDB UniProt sequences", unit="protein"):
        for header, seq in entries:
            # Extract protein ID from different header formats
            if "|" in header:
                # UniProt format: sp|PROTEINID|... or tr|PROTEINID|...
                extracted_id = header.split("|")[1]
            else:
                # Simple format: PROTEINID other_info
                extracted_id = header.split()[0]
            
            if protein_id == extracted_id:
                rows.append({
                    "Protein_ID": protein_id,
                    "IEDB_Sequence": seq
                })
                break
    
    return pd.DataFrame(rows)


matched_iedb_ids = (
    perfect_match["IEDB_Protein_ID"]
    .dropna()
    .unique()
    .tolist()
)

print("Fetching UniProt sequences for matched IEDB proteins...")
uniprot_seq_df = fetch_iedb_sequences(matched_iedb_ids)
print(f"Fetched sequences for {len(uniprot_seq_df)} IEDB proteins.")

# ---------------------- Prepare merge tables ---------------------- #
# Pathogen metadata: keep one row per pathogen protein ID
pathogen_meta = (
    pathogen_data[["Protein_ID", "Genus_species", "Sequence", "Strain"]]
    .drop_duplicates(subset=["Protein_ID"])
    .rename(columns={
        "Protein_ID": "Pathogen_Protein_ID_merge",
        "Genus_species": "Organism_Source",
        "Sequence": "Pathogen_Sequence"
    })
)

# IEDB source protein sequences: keep one row per IEDB protein ID
iedb_seq = (
    uniprot_seq_df
    .drop_duplicates(subset=["Protein_ID"])
    .rename(columns={
        "Protein_ID": "IEDB_Protein_ID_merge"
    })
)

print(f"Unique IEDB proteins with sequences: {len(iedb_seq)}")

# Keep one row per unique protein-protein pair
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

print(f"Unique protein-protein pairs for alignment: {len(full_align_analysis)}")

# Print unique pathogen protein IDs in the final dataset
unique_pathogen_proteins = full_align_analysis["Pathogen_Protein_ID"].nunique

# Alignment setup
aligner = Align.PairwiseAligner()
aligner.mode = "local"

# minimal penalties to avoid huge numbers of equivalent alignments
aligner.mismatch_score = -1
aligner.open_gap_score = -1
aligner.extend_gap_score = -0.5

def compute_identity_and_coverage(seq1, seq2):
    if pd.isna(seq1) or pd.isna(seq2):
        return pd.Series({
            "Percent_Identity": None,
            "Coverage": None
        })

    seq1 = str(seq1).strip()
    seq2 = str(seq2).strip()

    if not seq1 or not seq2:
        return pd.Series({
            "Percent_Identity": None,
            "Coverage": None
        })

    try:
        best = aligner.align(seq1, seq2)[0]
    except Exception:
        return pd.Series({
            "Percent_Identity": None,
            "Coverage": None
        })

    matches = 0
    aligned_length = 0

    for (start1, end1), (start2, end2) in zip(best.aligned[0], best.aligned[1]):
        segment1 = seq1[start1:end1]
        segment2 = seq2[start2:end2]

        for a, b in zip(segment1, segment2):
            aligned_length += 1
            if a == b:
                matches += 1

    if aligned_length == 0:
        return pd.Series({
            "Percent_Identity": 0.0,
            "Coverage": 0.0
        })

    percent_identity = (matches / aligned_length) * 100
    coverage = (aligned_length / min(len(seq1), len(seq2))) * 100

    return pd.Series({
        "Percent_Identity": percent_identity,
        "Coverage": coverage
    })

print("Computing percent identity and coverage...")

metrics = full_align_analysis.progress_apply(
    lambda row: compute_identity_and_coverage(
        row["IEDB_Sequence"],
        row["Pathogen_Sequence"]
    ),
    axis=1
)

full_align_analysis = pd.concat([full_align_analysis, metrics], axis=1)

# ---------------------- Summary table ---------------------- #
#%%

summary_df = (
    full_align_analysis
    .dropna(subset=["Percent_Identity", "Coverage"])
    .groupby(["IEDB_Protein_ID", "Epitope_Source"])
    .agg(
        Mean_Identity=("Percent_Identity", "mean"),
        SD_Identity=("Percent_Identity", "std"),
        Mean_Coverage=("Coverage", "mean"),
        SD_Coverage=("Coverage", "std"),
        N_matches=("IEDB_Protein_ID", "count")
    )
    .reset_index()
)

# replace NaN SD values for proteins with only one match
summary_df["SD_Identity"] = summary_df["SD_Identity"].fillna(0)
summary_df["SD_Coverage"] = summary_df["SD_Coverage"].fillna(0)

# optional: sort by number of matches
summary_df = summary_df.sort_values("N_matches", ascending=False)

summary_df = summary_df.round(2)

save_csv(summary_df, DATA_DIR / "proccesed/iedb_protein_similarity_summary.csv")

# ---------------------- Scatter plot ---------------------- #
#%%
plt.figure(figsize=(9, 7))

sc = plt.scatter(
    summary_df["Mean_Coverage"],
    summary_df["Mean_Identity"],
    c=summary_df["N_matches"],
    norm=LogNorm(),
    alpha=0.85
)

# error bars for SD
plt.errorbar(
    summary_df["Mean_Coverage"],
    summary_df["Mean_Identity"],
    xerr=summary_df["SD_Coverage"],
    yerr=summary_df["SD_Identity"],
    fmt="none",
    alpha=0.25,
    capsize=2
)

highlight_ids = [
    "P10809",
    "P0DMV8",
    "P11021",
    "P01308",
    "P62805",
    "P38646",
    "Q71DI3",
    "P62807",
    "P11021",
    "P11142",
    "P55087",
    "Q05329",
    "P06733"
]

label_df = summary_df[summary_df["IEDB_Protein_ID"].isin(highlight_ids)]

texts = []
for _, row in label_df.iterrows():
    txt = plt.text(
        row["Mean_Coverage"],
        row["Mean_Identity"],
        row["Epitope_Source"],
        fontsize=8
    )
    texts.append(txt)

adjust_text(texts, arrowprops=dict(arrowstyle='->', color='gray', lw=0.5))

plt.xlabel("Mean coverage (% of shorter protein)")
plt.ylabel("Mean percent identity (%)")

# colorbar
cbar = plt.colorbar(sc)
cbar.set_label("Number of matched pathogen proteins (log scale)")
plt.tight_layout()
plt.show()
# %%
