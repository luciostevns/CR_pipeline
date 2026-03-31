print("Starting script...")

import pandas as pd
from helpers import load_tsv, save_csv, load_tsv_robust, DATA_DIR

# Pathogen reference data
print("Loading isolates...")
pathogen_ref = load_tsv_robust(DATA_DIR / "raw/isolates.tsv")

# Proteome data
print("Loading UniProt proteomes...")
referenced_proteome = load_tsv(DATA_DIR / "raw/Uniprot_raw_refferenced_bacterial_proteomes.tsv")
other_proteome = load_tsv(DATA_DIR / "raw/Uniprot_raw_unrefferenced_bacterial_proteomes.tsv")

# Bind rows
proteome = pd.concat([referenced_proteome, other_proteome], ignore_index=True)

# Remove missing IDs needed for merge
proteome = proteome.dropna(subset=["Proteome Id", "Genome assembly ID"]).copy()
pathogen_ref = pathogen_ref.dropna(subset=["Assembly"]).copy()

# Merge proteome data onto pathogen reference by assembly
merged = proteome.merge(
    pathogen_ref,
    how="left",
    left_on="Genome assembly ID",
    right_on="Assembly",
    suffixes=("", "_ref")
)

# Keep only proteomes that matched a pathogenic isolate assembly
unique_merged = merged.dropna(subset=["Assembly"]).copy()

# Remove duplicate proteome IDs if present
unique_merged = unique_merged.drop_duplicates(subset=["Proteome Id"]).copy()

print(f"Matched pathogenic proteomes: {len(unique_merged)}")

# Save a clean proteome metadata table based on UniProt fields
proteome_metadata = unique_merged[
    ["Proteome Id", "Organism Id", "Genome assembly ID", "Protein count"]
].copy()

proteome_metadata = proteome_metadata.rename(columns={
    "Proteome Id": "Proteome_ID",
    "Organism Id": "Organism_ID",
    "Genome assembly ID": "Genome_Assembly_ID",
    "Protein count": "Protein_Count"
})

save_csv(
    proteome_metadata,
    DATA_DIR / "intermediate/pathogenic_bacteria_proteome.csv"
)

# Save proteome IDs for download of full proteomes
save_csv(
    proteome_metadata[["Proteome_ID"]],
    DATA_DIR / "intermediate/proteome_ids.csv"
)

# Diagnostics
print("Total annotated protein count:",
      pd.to_numeric(proteome_metadata["Protein_Count"], errors="coerce").sum())

matching_ids = set(proteome["Genome assembly ID"]) & set(pathogen_ref["Assembly"])
print(f"Number of matching Assembly IDs: {len(matching_ids)}")

print("Saved:")
print("- intermediate/pathogenic_bacteria_proteome.csv")
print("- intermediate/proteome_ids.csv")