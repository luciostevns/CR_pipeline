#%% IMPORTS
import pandas as pd
import requests
import time
import re

#%% LOAD DATA
perfect_match = pd.read_csv("../Data/perfect_matches_2_0.csv")
IEDB_data = pd.read_csv("../Data/wrangled_IEDB_with_sequences.csv")
cell_location_ref = pd.read_csv("../Data/Uniprot_subcellular_location_ref.csv")

#%% LOAD MULTIPLE DEEPLOCPRO NEGATIVE FILES
deeploc_neg1 = pd.read_csv("../Data/deeplocpro_Negative_1.csv")
deeploc_neg2 = pd.read_csv("../Data/deeplocpro_Negative_2.csv")
deeploc_neg3 = pd.read_csv("../Data/deeplocpro_Negative_3.csv")
deeploc_neg4 = pd.read_csv("../Data/deeplocpro_Negative_4.csv")
deeplocpro_positive = pd.read_csv("../Data/deeplocpro_Positive.csv")
deeplocpro_positive_2 = pd.read_csv("../Data/deeplocpro_Positive_2.csv")

deeplocpro_combined = pd.concat([deeploc_neg1, deeploc_neg2, deeploc_neg3, deeploc_neg4, deeplocpro_positive, deeplocpro_positive_2], ignore_index=True)

#%% CLEAN DEEPLOC PATHOGEN LOCATIONS (no mapping)
deeplocpro_combined.columns = deeplocpro_combined.columns.str.strip()
deeplocpro_combined["Localization"] = deeplocpro_combined["Localization"].str.strip()

deeploc_clean = deeplocpro_combined[["ACC", "Localization"]].rename(
    columns={"ACC": "Entry", "Localization": "pathogen_deeploc_subcellular_location"}
)

#%% CLEAN UNIPROT PATHOGEN LOCATIONS (no mapping)
uniprot_clean = cell_location_ref[["Entry", "Subcellular location [CC]"]].copy()
uniprot_clean["pathogen_uniprot_subcellular_location"] = (
    uniprot_clean["Subcellular location [CC]"]
    .str.extract(r"SUBCELLULAR LOCATION:\s*([^;{\.]+)", flags=re.IGNORECASE)[0]
    .str.strip()
)
uniprot_clean = uniprot_clean.drop(columns=["Subcellular location [CC]"])

#%% MERGE PATHOGEN LOCATION INFO INTO perfect_match
perfect_match["Pathogen_Protein_ID"] = perfect_match["Pathogen_Protein_ID"].str.strip().str.upper()
uniprot_clean["Entry"] = uniprot_clean["Entry"].str.strip().str.upper()
deeploc_clean["Entry"] = deeploc_clean["Entry"].str.strip().str.upper()

perfect_match = perfect_match.merge(uniprot_clean, how="left", left_on="Pathogen_Protein_ID", right_on="Entry").drop(columns=["Entry"])
perfect_match = perfect_match.merge(deeploc_clean, how="left", left_on="Pathogen_Protein_ID", right_on="Entry").drop(columns=["Entry"])

#%% FETCH EPIOTOPE UNIPROT LOCATIONS VIA API (no mapping)
epitope_ids = IEDB_data["Protein_ID"].dropna().unique()

def fetch_uniprot_location(uniprot_id):
    """
    Search UniProt REST API for subcellular location of a given UniProt ID.
    Returns the subcellular location as a string, or None if not found.
    """

    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    try:
        response = requests.get(url)
        if response.status_code != 200:
            return None
        data = response.json()
        for comment in data.get("comments", []):
            if comment.get("commentType") == "SUBCELLULAR LOCATION":
                locations = [
                    loc.get("location", {}).get("value")
                    for loc in comment.get("subcellularLocations", [])
                    if loc.get("location")
                ]
                return "; ".join(locations) if locations else None
        return None
    except Exception as e:
        print(f"Error fetching {uniprot_id}: {e}")
        return None

epitope_location_map = {}
for uid in epitope_ids:
    epitope_location_map[uid] = fetch_uniprot_location(uid)
    time.sleep(0.5)

# Apply to both datasets
IEDB_data["epitope_uniprot_subcellular_location"] = IEDB_data["Protein_ID"].map(epitope_location_map)
perfect_match["epitope_uniprot_subcellular_location"] = perfect_match["IEDB_Protein_ID"].map(epitope_location_map)

#%% MERGE DEEPLOC EPIOTOPE LOCATION (from 2 files, no mapping)
epitope_deeploc1 = pd.read_csv("../Data/deeploc_epitopes_1.csv")
epitope_deeploc2 = pd.read_csv("../Data/deeploc_epitopes_2.csv")

epitope_deeploc = pd.concat([epitope_deeploc1, epitope_deeploc2], ignore_index=True)
epitope_deeploc.columns = epitope_deeploc.columns.str.strip()
epitope_deeploc["Localizations"] = epitope_deeploc["Localizations"].str.strip()

epitope_deeploc_clean = epitope_deeploc[["Protein_ID", "Localizations"]].rename(
    columns={"Localizations": "epitope_deeploc_subcellular_location"}
)

# Merge into both datasets
IEDB_data = IEDB_data.merge(
    epitope_deeploc_clean,
    how="left",
    on="Protein_ID"
)

perfect_match = perfect_match.merge(
    epitope_deeploc_clean,
    how="left",
    left_on="IEDB_Protein_ID",
    right_on="Protein_ID"
).drop(columns=["Protein_ID"])

#%% SAVE
perfect_match.to_csv("../Data/perfect_matches_finished.csv", index=False)
IEDB_data.to_csv("../Data/IEDB_with_locations.csv", index=False)
print("✅ Saved finished files with raw UniProt and DeepLoc locations to:")
print("   → perfect_matches_finished.csv")
print("   → IEDB_with_locations.csv")
