#%%
import pandas as pd
import requests
from Bio import Align
from tqdm import tqdm
import io

# ---------------------- Step 1: Load Data ---------------------- #
perfect_match = pd.read_csv("../Data/perfect_matches_finished.csv")
pathogen_data = pd.read_csv("../Data/wrangled_all_pathogen_prots.csv")

# Filter out unmatched rows (no matched 9mer = no match)
perfect_match = perfect_match[perfect_match["Matched_9mer"].notna()].copy()

# Normalize IDs
perfect_match["Pathogen_Protein_ID"] = perfect_match["Pathogen_Protein_ID"].astype(str).str.strip().str.upper()
perfect_match["IEDB_Protein_ID"] = perfect_match["IEDB_Protein_ID"].astype(str).str.strip().str.upper()
pathogen_data["Protein_ID"] = pathogen_data["Protein_ID"].astype(str).str.strip().str.upper()

# Get unique IEDB protein IDs that were matched
matched_iedb_ids = perfect_match["IEDB_Protein_ID"].dropna().unique().tolist()

# ---------------------- Step 2: Fetch UniProt Sequences ---------------------- #
def fetch_uniprot_sequences(uniprot_ids):
    url = "https://rest.uniprot.org/uniprotkb/stream"
    headers = {"Accept": "text/tab-separated-values"}
    query = f"({' OR '.join(uniprot_ids)})"
    params = {
        "query": f"accession:({query})",
        "fields": "accession,sequence",
        "format": "tsv"
    }

    response = requests.get(url, headers=headers, params=params)
    response.raise_for_status()

    df = pd.read_csv(io.StringIO(response.text), sep="\t")
    return df.rename(columns={"Entry": "Protein_ID", "Sequence": "Full_Sequence"})

print("Fetching UniProt sequences...")
uniprot_seq_df = fetch_uniprot_sequences(matched_iedb_ids)
print(f"Fetched sequences for {len(uniprot_seq_df)} IEDB proteins.")

# ---------------------- Step 3: Rename and Merge ---------------------- #
# Rename to avoid conflict with 'Sequence' column in perfect_match
pathogen_data = pathogen_data.rename(columns={"Sequence": "Pathogen_Sequence"})

# Merge sequences
full_align_analysis = (
    perfect_match
    .drop_duplicates(subset=["IEDB_Protein_ID", "Pathogen_Protein_ID"])
    .merge(
        pathogen_data[["Protein_ID", "Genus_Species", "Pathogen_Sequence"]],
        how="left",
        left_on="Pathogen_Protein_ID",
        right_on="Protein_ID"
    )
    .merge(
        uniprot_seq_df.drop_duplicates("Protein_ID"),
        how="left",
        left_on="IEDB_Protein_ID",
        right_on="Protein_ID"
    )
)

# Rename for consistency
full_align_analysis = full_align_analysis.rename(columns={
    "Full_Sequence": "IEDB_Sequence",
    "Genus_Species": "Organism_Source"
})

# ---------------------- Step 4: Compute Alignment ---------------------- #
aligner = Align.PairwiseAligner()
aligner.mode = 'global'

def compute_similarity(seq1, seq2):
    if pd.isna(seq1) or pd.isna(seq2):
        return None
    score = aligner.score(seq1, seq2)
    max_length = max(len(seq1), len(seq2))
    return (score / max_length) * 100 if max_length else 0

percent_identities = []
for idx, row in tqdm(full_align_analysis.iterrows(), total=len(full_align_analysis), desc="Computing Percent Identity"):
    pid = compute_similarity(row["IEDB_Sequence"], row["Pathogen_Sequence"])
    percent_identities.append(pid)

full_align_analysis["Percent_Identity"] = percent_identities

# ---------------------- Step 5: Export ---------------------- #
result_with_similarity = full_align_analysis[[
    "IEDB_Protein_ID",
    "Pathogen_Protein_ID",
    "Strain",
    "Organism_Source",
    "Assay_ID",
    "Epitope_Source",
    "Disease",
    "Percent_Identity"
]]

result_with_similarity.to_csv("../Data/full_align_with_similarity.csv", index=False)
print("\nFinished! Exported to ../Data/full_align_with_similarity.csv ✅")

#%%
