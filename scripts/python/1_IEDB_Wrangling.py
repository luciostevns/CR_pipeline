# Imports
import pandas as pd
import re
from collections import defaultdict
from helpers import load_excel, save_csv, DATA_DIR

################ Initial filtering and wrangling ###############################

# Read the full IEDB dataset
autoimmune_data = load_excel(DATA_DIR / "raw/IEDB_autoimmune_epitope_assays_raw.xlsx")
print(f"Original dataset size: {len(autoimmune_data)}")

# Rename columns
autoimmune_data = autoimmune_data.rename(columns={
    "Assay ID - IEDB IRI": "Assay_ID",
    "Epitope - Name": "Sequence",
    "Epitope - Molecule Parent": "Protein_source",
    "Epitope - Molecule Parent IRI": "Protein_ID",
    "1st in vivo Process - Disease": "Disease",
    "1st in vivo Process - Disease Stage": "Disease_stage",
    "MHC Restriction - Name": "MHC_restriction",
    "Epitope - Starting Position": "epitope_start_pos",
    "Epitope - Ending Position": "epitope_end_pos"
})

# Clean sequence and protein ID
autoimmune_data["Sequence"] = (
    autoimmune_data["Sequence"]
    .str.replace(r"\+.*", "", regex=True)
    .str.strip()
)

autoimmune_data["Protein_ID"] = autoimmune_data["Protein_ID"].str.extract(r"([^/]+$)")

# Filtering
autoimmune_data_wrangled = autoimmune_data[
    autoimmune_data["Epitope - Modified residues"].isna() &
    autoimmune_data["Sequence"].str.len().between(12, 25)
].copy()
print(f"Dataset size after basic filtering: {len(autoimmune_data_wrangled)}")

# Remove duplicate sequences with same protein ID (keep first occurrence)
autoimmune_data_wrangled = autoimmune_data_wrangled.drop_duplicates(
    subset=["Protein_ID", "Sequence"]
)
print(f"Dataset size after removing duplicate seqs: {len(autoimmune_data_wrangled)}")

# Select relevant columns
autoimmune_data_wrangled = autoimmune_data_wrangled[
    [
        "Assay_ID",
        "Sequence",
        "Protein_ID",
        "Protein_source",
        "Disease",
        "Disease_stage",
        "MHC_restriction",
        "epitope_start_pos",
        "epitope_end_pos"
    ]
]

################ Nested epitope removal #######################################

grouped_results = []

for protein_id, group in autoimmune_data_wrangled.groupby("Protein_ID", sort=False):
    group = group.copy()
    seqs = group["Sequence"].tolist()

    is_nested = []
    for seq in seqs:
        nested = any(
            (len(other) > len(seq)) and (seq in other)
            for other in seqs
        )
        is_nested.append(nested)

    group["Is_nested"] = is_nested
    group["Protein_ID"] = protein_id
    grouped_results.append(group)

autoimmune_data_wrangled = pd.concat(grouped_results, ignore_index=True)

filtered_sequences = autoimmune_data_wrangled[
    ~autoimmune_data_wrangled["Is_nested"]
].drop(columns=["Is_nested"])

print(f"Dataset size after removing nested epitopes: {len(filtered_sequences)}")
# Write output
print(filtered_sequences.columns.tolist())
save_csv(filtered_sequences, DATA_DIR / "intermediate/wrangled_IEDB.csv")