# Imports
import pandas as pd
import re
from collections import defaultdict
from helpers import load_excel, save_csv, DATA_DIR

################ Initial filtering and wrangling ###############################

# Read the full IEDB dataset
autoimmune_data = load_excel(DATA_DIR / "raw/IEDB_autoimmune_epitope_assays_raw.xlsx")

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
]

# Remove duplicate sequences (keep first occurrence)
autoimmune_data_wrangled = autoimmune_data_wrangled.drop_duplicates(
    subset="Sequence"
)

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

################ Nested proteins removal ######################################

# Function to generate all possible 9-mers from a sequence
def generate_9mers(seq):
    if len(seq) < 9:
        return []
    return [seq[i:i+9] for i in range(len(seq) - 8)]

# Generate 9-mers for each sequence
autoimmune_data_wrangled["Nine_mers"] = autoimmune_data_wrangled["Sequence"].apply(generate_9mers)

# Expand dataframe: one row per 9-mer
expanded_9mers = (
    autoimmune_data_wrangled[["Sequence", "Nine_mers"]]
    .explode("Nine_mers")
    .dropna()
)

# Precompute for speed
nine_mer_dict = defaultdict(list)

for nine_mer, seq_len in zip(
    expanded_9mers["Nine_mers"],
    expanded_9mers["Sequence"].str.len()
):
    nine_mer_dict[nine_mer].append(seq_len)

# Identify nested sequences
def is_nested(row):
    seq_len = len(row["Sequence"])

    for nine_mer in row["Nine_mers"]:
        if nine_mer in nine_mer_dict:
            for other_len in nine_mer_dict[nine_mer]:
                if other_len > seq_len:
                    return True

    return False

autoimmune_data_wrangled["Is_nested"] = autoimmune_data_wrangled.apply(
    is_nested,
    axis=1
)

# Keep only non-nested sequences
filtered_sequences = autoimmune_data_wrangled[
    ~autoimmune_data_wrangled["Is_nested"]
].drop(columns=["Nine_mers", "Is_nested"])

# Write output
save_csv(filtered_sequences, DATA_DIR / "intermediate/wrangled_IEDB.csv")
