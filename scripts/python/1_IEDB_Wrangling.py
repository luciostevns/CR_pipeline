#%%
import pandas as pd
import re
from helpers import load_tsv

################ Initial filtering and wrangling ###############################

# Read the full IEDB dataset
autoimmune_data = pd.read_csv(
    "./../Data/tcell_table_export_1740575887.tsv",
    sep="\t"
)

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
expanded_9mers_list = list(
    zip(expanded_9mers["Nine_mers"], expanded_9mers["Sequence"].str.len())
)

# Identify nested sequences
def is_nested(row):
    seq_len = len(row["Sequence"])
    for nine_mer in row["Nine_mers"]:
        for other_nine_mer, other_len in expanded_9mers_list:
            if nine_mer == other_nine_mer and other_len > seq_len:
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
filtered_sequences.to_csv(
    "./../Data/wrangled_IEDB.csv",
    index=False
)

# %%
