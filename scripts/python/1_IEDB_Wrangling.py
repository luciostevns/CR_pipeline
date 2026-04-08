# Imports
import pandas as pd
from helpers import load_excel, save_csv, DATA_DIR

################ Load datasets ###############################################

autoimmune_ref = load_excel(DATA_DIR / "raw/IEDB_autoimmune_epitope_assays_raw.xlsx")
general_IEDB_data = load_excel(DATA_DIR / "raw/IEDB_epitop_assays_raw.xlsx")

print(f"Autoimmune reference dataset size: {len(autoimmune_ref)}")
print(f"General IEDB dataset size: {len(general_IEDB_data)}")

################ Rename columns ##############################################

rename_dict = {
    "Assay ID - IEDB IRI": "Assay_ID",
    "Epitope - Name": "Sequence",
    "Epitope - Molecule Parent": "Protein_source",
    "Epitope - Molecule Parent IRI": "Protein_ID",
    "1st in vivo Process - Disease": "Disease",
    "Epitope - Starting Position": "epitope_start_pos",
    "Epitope - Ending Position": "epitope_end_pos",
    "Epitope - Modified residues": "Modified_residues"
}

cols_to_keep = list(rename_dict.values())

autoimmune_ref = (
    autoimmune_ref
    .rename(columns=rename_dict)
    .filter(cols_to_keep)
)

general_IEDB_data = (
    general_IEDB_data
    .rename(columns=rename_dict)
    .filter(cols_to_keep)
)

################ Clean disease column ########################################

autoimmune_ref["Disease"] = (
    autoimmune_ref["Disease"]
    .fillna("")
    .astype(str)
    .str.strip()
)

general_IEDB_data["Disease"] = (
    general_IEDB_data["Disease"]
    .fillna("")
    .astype(str)
    .str.strip()
)

################ Add missing autoimmune diseases #############################

extra_autoimmune_diseases = {
    "type 1 diabetes mellitus"
}

extra_data = general_IEDB_data[
    general_IEDB_data["Disease"].isin(extra_autoimmune_diseases)
].copy()

print(f"Extra autoimmune rows added: {len(extra_data)}")
print(extra_data["Disease"].value_counts())

################ Combine datasets ############################################

combined_data = pd.concat([autoimmune_ref, extra_data], ignore_index=True)

# Remove duplicates (important!)
combined_data = combined_data.drop_duplicates(subset=["Assay_ID"])

print(f"Combined dataset size: {len(combined_data)}")

################ Basic cleaning ##############################################

combined_data["Sequence"] = (
    combined_data["Sequence"]
    .fillna("")
    .astype(str)
    .str.replace(r"\+.*", "", regex=True)
    .str.strip()
)

combined_data["Protein_ID"] = (
    combined_data["Protein_ID"]
    .fillna("")
    .astype(str)
    .str.extract(r"([^/]+$)", expand=False)
)

################ Filtering ###################################################

combined_data_wrangled = combined_data[
    combined_data["Modified_residues"].isna() &
    combined_data["Sequence"].str.len().between(12, 25, inclusive="both")
].copy()

print(f"Dataset size after basic filtering: {len(combined_data_wrangled)}")

# Remove duplicate sequences per protein
combined_data_wrangled = combined_data_wrangled.drop_duplicates(
    subset=["Protein_ID", "Sequence"]
)

print(f"Dataset size after removing duplicate seqs: {len(combined_data_wrangled)}")

################ Nested epitope removal ######################################

grouped_results = []

for protein_id, group in combined_data_wrangled.groupby("Protein_ID", sort=False):
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
    grouped_results.append(group)

combined_data_wrangled = pd.concat(grouped_results, ignore_index=True)

filtered_sequences = combined_data_wrangled[
    ~combined_data_wrangled["Is_nested"]
].drop(columns=["Is_nested"])

print(f"Dataset size after removing nested epitopes: {len(filtered_sequences)}")

################ Save #########################################################
save_csv(filtered_sequences, DATA_DIR / "intermediate/wrangled_IEDB.csv")

