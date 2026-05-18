import pandas as pd
import numpy as np
from pathlib import Path
from helpers import load_excel, save_csv, load_csv, DATA_DIR

# -----------------------------
# Load perfect match results
# -----------------------------
perfect_match = load_csv(DATA_DIR / "proccesed/perfect_matches_2_0.csv")

perfect_match["Matched"] = (
    perfect_match["Matched"]
    .astype(str)
    .str.lower()
    .eq("true")
)

# -----------------------------
# NetSurfP column names
# -----------------------------
netsurfp_cols = [
    "id", "seq", "n", "rsa", "asa", "q3",
    "p[q3_H]", "p[q3_E]", "p[q3_C]",
    "q8", "p[q8_G]", "p[q8_H]", "p[q8_I]", "p[q8_B]",
    "p[q8_E]", "p[q8_S]", "p[q8_T]", "p[q8_C]",
    "phi", "psi", "disorder"
]


# -----------------------------
# Read NetSurfP IEDB manually
# -----------------------------
netsurfp_dir = DATA_DIR / "raw" / "netsurfp"

rows = []

with open(netsurfp_dir / "netsurfp_IEDB_prediction.csv", encoding="utf-8") as f:
    next(f)  # skip header

    for line_number, line in enumerate(f, start=2):
        line = line.strip()

        if not line:
            continue

        parts = line.rsplit(",", 20)

        if len(parts) != 21:
            print(f"Skipping malformed line {line_number}: found {len(parts)} fields")
            continue

        rows.append(parts)

netsurfp_iedb = pd.DataFrame(rows, columns=netsurfp_cols)


# -----------------------------
# Read NetSurfP pathogen files manually
# -----------------------------
netsurfp_pathogen_files = sorted(
    netsurfp_dir.glob("netsurfp_pathogen_prediction_*.csv")
)

pathogen_rows = []

for file in netsurfp_pathogen_files:
    with open(file, encoding="utf-8") as f:
        next(f)  # skip header

        for line_number, line in enumerate(f, start=2):
            line = line.strip()

            if not line:
                continue

            parts = line.rsplit(",", 20)

            if len(parts) != 21:
                print(f"Skipping malformed line {file.name}, line {line_number}: found {len(parts)} fields")
                continue

            pathogen_rows.append(parts)

netsurfp_pathogen = pd.DataFrame(pathogen_rows, columns=netsurfp_cols)

# -----------------------------
# Read DeepLoc files normally
# -----------------------------
deeploc_dir = DATA_DIR / "raw" / "deeploc"

deeploc_iedb = pd.read_csv(
    deeploc_dir / "deeploc_IEDB_prediction.csv",
    dtype=str
)

deeploc_iedb.columns = deeploc_iedb.columns.str.strip()

deeploc_pathogen_files = sorted(
    deeploc_dir.glob("deeplocpro_prediction_*.csv")
)

deeploc_pathogen = pd.concat(
    [
        pd.read_csv(f, dtype=str)
        for f in deeploc_pathogen_files
    ],
    ignore_index=True
)

deeploc_pathogen.columns = deeploc_pathogen.columns.str.strip()

# Clean protein ID vals for check
# perfect_match
perfect_match["Pathogen_Protein_ID_clean"] = (
    perfect_match["Pathogen_Protein_ID"]
    .astype("string")
    .str.strip()
    .str.upper()
)

perfect_match["IEDB_Protein_ID_clean"] = (
    perfect_match["IEDB_Protein_ID"]
    .astype("string")
    .str.strip()
    .str.upper()
)

# NetSurfP
netsurfp_iedb["Protein_ID_clean"] = (
    netsurfp_iedb["id"]
    .astype("string")
    .str.extract(r">(?:sp|tr)_([^_]+)_", expand=False)
    .str.upper()
)

netsurfp_pathogen["Protein_ID_clean"] = (
    netsurfp_pathogen["id"]
    .astype(str)
    .str.extract(r"^>?([A-Za-z0-9]+)", expand=False)
    .str.upper()
)

# DeepLoc
deeploc_iedb["Protein_ID_clean"] = (
    deeploc_iedb["Protein_ID"]
    .astype(str)
    .str.extract(r"(?:sp|tr)_([^_]+)_", expand=False)
    .str.upper()
)

deeploc_pathogen["Protein_ID_clean"] = (
    deeploc_pathogen["ACC"]
    .astype(str)
    .str.strip()
    .str.upper()
)

# Remove index-like columns if present
deeploc_pathogen = deeploc_pathogen.loc[:, ~deeploc_pathogen.columns.str.contains("^Unnamed")]
deeploc_iedb = deeploc_iedb.loc[:, ~deeploc_iedb.columns.str.contains("^Unnamed")]
netsurfp_pathogen = netsurfp_pathogen.loc[:, ~netsurfp_pathogen.columns.str.contains("^Unnamed")]
netsurfp_iedb = netsurfp_iedb.loc[:, ~netsurfp_iedb.columns.str.contains("^Unnamed")]

# Check for proteins ID
# -----------------------------
# Expected protein IDs from matched rows
# -----------------------------
matched_iedb_ids = set(
    perfect_match.loc[perfect_match["Matched"], "IEDB_Protein_ID_clean"].dropna()
)

matched_pathogen_ids = set(
    perfect_match.loc[perfect_match["Matched"], "Pathogen_Protein_ID_clean"].dropna()
)

# -----------------------------
# Observed protein IDs in predictions
# -----------------------------
netsurfp_iedb_ids = set(netsurfp_iedb["Protein_ID_clean"].dropna())
netsurfp_pathogen_ids = set(netsurfp_pathogen["Protein_ID_clean"].dropna())

deeploc_iedb_ids = set(deeploc_iedb["Protein_ID_clean"].dropna())
deeploc_pathogen_ids = set(deeploc_pathogen["Protein_ID_clean"].dropna())

# -----------------------------
# Missing checks
# -----------------------------
missing_netsurfp_iedb = matched_iedb_ids - netsurfp_iedb_ids
missing_netsurfp_pathogen = matched_pathogen_ids - netsurfp_pathogen_ids

missing_deeploc_iedb = matched_iedb_ids - deeploc_iedb_ids
missing_deeploc_pathogen = matched_pathogen_ids - deeploc_pathogen_ids

print("Expected matched IEDB proteins:", len(matched_iedb_ids))
print("Expected matched pathogen proteins:", len(matched_pathogen_ids))

print("\nMissing from NetSurfP IEDB:", len(missing_netsurfp_iedb))
print("Missing from NetSurfP pathogen:", len(missing_netsurfp_pathogen))

print("\nMissing from DeepLoc IEDB:", len(missing_deeploc_iedb))
print("Missing from DeepLoc pathogen:", len(missing_deeploc_pathogen))

# Extract missing ids
# Convert missing pathogen NetSurfP IDs to a sorted dataframe
missing_netsurfp_pathogen_df = pd.DataFrame({
    "Pathogen_Protein_ID": sorted(missing_netsurfp_pathogen)
})

print(missing_netsurfp_pathogen_df.shape)
print(missing_netsurfp_pathogen_df.head())

# Save missing IDs
save_csv(
    missing_netsurfp_pathogen_df,
    DATA_DIR / "proccesed/missing_netsurfp_pathogen_ids.csv"
)

missing_netsurfp_pathogen_meta = (
    perfect_match[
        perfect_match["Matched"] &
        perfect_match["Pathogen_Protein_ID_clean"].isin(missing_netsurfp_pathogen)
    ][
        [
            "Pathogen_Protein_ID",
            "Pathogen_Organism",
            "Pathogen_Scientific_name",
            "Pathogen_Strain",
            "Pathogen_Annotation",
            "Pathogen_Gene_Name"
        ]
    ]
    .drop_duplicates()
    .sort_values("Pathogen_Protein_ID")
)

print(missing_netsurfp_pathogen_meta.shape)
print(missing_netsurfp_pathogen_meta.head())

# Add deeploc and netsurfP meta data to match df
# ---- DeepLoc metadata ----
iedb_location_meta = (
    deeploc_iedb[["Protein_ID_clean", "Localizations"]]
    .dropna(subset=["Protein_ID_clean"])
    .drop_duplicates("Protein_ID_clean")
    .rename(columns={"Localizations": "IEDB_prot_location"})
)

pathogen_location_meta = (
    deeploc_pathogen[["Protein_ID_clean", "Localization"]]
    .dropna(subset=["Protein_ID_clean"])
    .drop_duplicates("Protein_ID_clean")
    .rename(columns={"Localization": "pathogen_prot_location"})
)

# ---- Merge DeepLoc metadata ----
perfect_match_meta = (
    perfect_match
    .merge(
        pathogen_location_meta,
        how="left",
        left_on="Pathogen_Protein_ID_clean",
        right_on="Protein_ID_clean"
    )
    .drop(columns=["Protein_ID_clean"])
    .merge(
        iedb_location_meta,
        how="left",
        left_on="IEDB_Protein_ID_clean",
        right_on="Protein_ID_clean"
    )
    .drop(columns=["Protein_ID_clean"])
)

# ---- Make numeric columns ----
for df in [netsurfp_iedb, netsurfp_pathogen]:
    df["n"] = pd.to_numeric(df["n"], errors="coerce")
    df["rsa"] = pd.to_numeric(df["rsa"], errors="coerce")

for col in ["Pathogen_Start", "Pathogen_End", "Epitope_Start", "Epitope_End"]:
    perfect_match_meta[col] = pd.to_numeric(perfect_match_meta[col], errors="coerce")

# ---- Give each match row a unique id ----
perfect_match_meta = perfect_match_meta.reset_index(drop=True)
perfect_match_meta["match_row_id"] = perfect_match_meta.index

# ---- Pathogen matched peptide RSA ----
pathogen_regions = perfect_match_meta[
    ["match_row_id", "Pathogen_Protein_ID_clean", "Pathogen_Start", "Pathogen_End"]
].dropna()

pathogen_rsa = pathogen_regions.merge(
    netsurfp_pathogen[["Protein_ID_clean", "n", "rsa"]],
    how="left",
    left_on="Pathogen_Protein_ID_clean",
    right_on="Protein_ID_clean"
)

pathogen_rsa = pathogen_rsa[
    pathogen_rsa["n"].between(
        pathogen_rsa["Pathogen_Start"],
        pathogen_rsa["Pathogen_End"]
    )
]

pathogen_rsa_summary = (
    pathogen_rsa
    .groupby("match_row_id")["rsa"]
    .mean()
    .reset_index(name="rsa_pathogen_peptide_mean")
)

# ---- IEDB epitope RSA ----
iedb_regions = perfect_match_meta[
    ["match_row_id", "IEDB_Protein_ID_clean", "Epitope_Start", "Epitope_End"]
].dropna()

iedb_rsa = iedb_regions.merge(
    netsurfp_iedb[["Protein_ID_clean", "n", "rsa"]],
    how="left",
    left_on="IEDB_Protein_ID_clean",
    right_on="Protein_ID_clean"
)

iedb_rsa = iedb_rsa[
    iedb_rsa["n"].between(
        iedb_rsa["Epitope_Start"],
        iedb_rsa["Epitope_End"]
    )
]

iedb_rsa_summary = (
    iedb_rsa
    .groupby("match_row_id")["rsa"]
    .mean()
    .reset_index(name="rsa_epitope_mean")
)

# ---- Merge RSA summaries back ----
perfect_match_meta = (
    perfect_match_meta
    .merge(pathogen_rsa_summary, how="left", on="match_row_id")
    .merge(iedb_rsa_summary, how="left", on="match_row_id")
).drop(columns=["match_row_id"])

# ---- Quick diagnostics ----
print("Rows:", len(perfect_match_meta))

print("\nDeepLoc merge coverage:")
print("Pathogen location missing:", perfect_match_meta["pathogen_prot_location"].isna().sum())
print("IEDB location missing:", perfect_match_meta["IEDB_prot_location"].isna().sum())

print("\nNetSurfP RSA coverage:")
print("Pathogen RSA missing:", perfect_match_meta["rsa_pathogen_peptide_mean"].isna().sum())
print("IEDB RSA missing:", perfect_match_meta["rsa_epitope_mean"].isna().sum())

# Save cleaned prediction outputs
save_csv(perfect_match_meta, DATA_DIR / "proccesed/perfect_match_DN_add.csv")
print(perfect_match_meta["Matched"].value_counts(dropna=False))