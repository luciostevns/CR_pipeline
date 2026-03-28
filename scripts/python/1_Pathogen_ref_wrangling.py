#%%
import pandas as pd
import numpy as np
    
# Pathogen reference data
pathogen_ref = pd.read_csv("../Data/pathogen_ref.tsv", sep='\t')

# Proteome data
referenced_proteome = pd.read_csv("../Data/referenced_proteomes.tsv", sep='\t')
other_proteome = pd.read_csv("../Data/other_proteomes.tsv", sep='\t')

# binding rows of ref and other
proteome = pd.concat([referenced_proteome, other_proteome])

# Remove NA's
proteome = proteome.dropna(subset=["Proteome Id"])
proteome = proteome.dropna(subset=["Genome assembly ID"])

# Pattern to dettect
pattern = r'(?i)human|homo'

# Remove NA's
pathogen_ref = pathogen_ref.dropna(subset=['Host'])
pathogen_ref = pathogen_ref.dropna(subset=['Assembly'])

# Filter out rows that don't contain 'human' in the 'Host' column
filtered_df = pathogen_ref[pathogen_ref['Host'].str.contains(pattern, regex=True)]

# Exclude rows that contain 'non-human' in the 'Host' column
filtered_pathogen_ref = filtered_df[~filtered_df['Host'].str.contains(r'(?i)non', regex=True)]

# Merging proteome data onto ref
merged_proteome = proteome.merge(filtered_pathogen_ref, how="left", left_on="Genome assembly ID", right_on="Assembly", suffixes=('', '_ref'))

# Remove weird instances with missing Proteime id
# Looks at NA's
unique_merged = merged_proteome.dropna(subset=['Assembly'])

# Saving certain columns from merged df as .csv
unique_merged[["Proteome Id", "#Organism group", "Strain", "Taxonomic lineage", "Protein count", "Assembly"]].to_csv("../Data/Pathogenic_bacteria_proteome.csv", index=False)


# Saving proteomes IDs as .txt for download of full proteomes
unique_merged['Proteome Id'].to_csv("../Data/proteome_ids.txt", header=False, index=False)

# Manually removing empty proteomes files
empty_proteomes = ["UP001298242","UP001299051","UP001299808","UP001304105", "UP001311976","UP001317873","UP001321007","UP001347461","UP001362035","UP001368540","UP001370664","UP001375347","UP001381400","UP001389104","UP001392488","UP001394715","UP001396551","UP001435425","UP001441349","UP001458050","UP001473423","UP001488889"]

# Finding annotated total protein count
print(unique_merged.shape)
print(unique_merged["Protein count"].sum())

# Check how many IDs match before merging
matching_ids = set(proteome["Genome assembly ID"]) & set(filtered_pathogen_ref["Assembly"])
print(f"Number of matching Assembly IDs: {len(matching_ids)}")





# %%
