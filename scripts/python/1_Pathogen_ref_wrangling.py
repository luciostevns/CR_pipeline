print("Starting script...")   
import pandas as pd
import numpy as np
from helpers import load_tsv, save_csv, load_tsv_robust, DATA_DIR

# Pathogen reference data
print("Loading isolates...")
pathogen_ref = load_tsv_robust(DATA_DIR / "raw/isolates.tsv")

# Proteome data
print("Loading Uniprot proteomes...")
referenced_proteome = load_tsv(DATA_DIR / "raw/Uniprot_raw_refferenced_bacterial_proteomes.tsv")
other_proteome = load_tsv(DATA_DIR / "raw/Uniprot_raw_unrefferenced_bacterial_proteomes.tsv")

# binding rows of ref and other
proteome = pd.concat([referenced_proteome, other_proteome])

# Remove NA's
proteome = proteome.dropna(subset=["Proteome Id"])
proteome = proteome.dropna(subset=["Genome assembly ID"])
pathogen_ref = pathogen_ref.dropna(subset=['Assembly'])

# Merging proteome data onto ref
merged = proteome.merge(pathogen_ref, how="left", left_on="Genome assembly ID", right_on="Assembly", suffixes=('', '_ref'))

# Remove NaN assembly ids from merged df
unique_merged = merged.dropna(subset=["Assembly"])

# Saving certain columns from merged df as .csv
save_csv(unique_merged[["Proteome Id", "#Organism group", "Strain", "TaxID", "Protein count", "Assembly"]], DATA_DIR / "intermediate/pathogenic_bacteria_proteome.csv")

# Saving proteomes IDs as .txt for download of full proteomes
save_csv(unique_merged["Proteome Id"], DATA_DIR / "intermediate/proteome_ids.csv")

# Finding annotated total protein count
print("protein count:", pd.to_numeric(unique_merged["Protein count"], errors="coerce").sum())

# Check how many IDs match before merging
matching_ids = set(proteome["Genome assembly ID"]) & set(pathogen_ref["Assembly"])
print(f"Number of matching Assembly IDs: {len(matching_ids)}")
