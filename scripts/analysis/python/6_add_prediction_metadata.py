import pandas as pd

from crossreactivity.io import DATA_DIR, load_csv, save_csv


NETSURFP_COLS = [
    "id", "seq", "n", "rsa", "asa", "q3",
    "p[q3_H]", "p[q3_E]", "p[q3_C]",
    "q8", "p[q8_G]", "p[q8_H]", "p[q8_I]", "p[q8_B]",
    "p[q8_E]", "p[q8_S]", "p[q8_T]", "p[q8_C]",
    "phi", "psi", "disorder",
]


def read_netsurfp_file(path) -> pd.DataFrame:
    rows = []

    with open(path, encoding="utf-8") as f:
        next(f)

        for line_number, line in enumerate(f, start=2):
            line = line.strip()

            if not line:
                continue

            parts = line.rsplit(",", 20)

            if len(parts) != 21:
                print(
                    f"Skipping malformed line {path.name}, "
                    f"line {line_number}: found {len(parts)} fields"
                )
                continue

            rows.append(parts)

    return pd.DataFrame(rows, columns=NETSURFP_COLS)


def load_prediction_inputs():
    netsurfp_dir = DATA_DIR / "raw" / "netsurfp"
    deeploc_dir = DATA_DIR / "raw" / "deeploc"

    print("Loading perfect match results...")
    perfect_match = load_csv(DATA_DIR / "proccesed/perfect_matches_2_0.csv")

    perfect_match["Matched"] = (
        perfect_match["Matched"]
        .astype(str)
        .str.lower()
        .eq("true")
    )

    print("Loading NetSurfP IEDB predictions...")
    netsurfp_iedb_files = sorted(
        netsurfp_dir.glob("netsurfp_IEDB_prediction_*.csv")
    )

    print("Loading NetSurfP pathogen predictions...")
    netsurfp_pathogen_files = sorted(
        netsurfp_dir.glob("netsurfp_pathogen_prediction_*.csv")
    )

    netsurfp_pathogen = pd.concat(
        [read_netsurfp_file(file) for file in netsurfp_pathogen_files],
        ignore_index=True
    )

    netsurfp_iedb = pd.concat(
        [read_netsurfp_file(file) for file in netsurfp_iedb_files],
        ignore_index=True
    )

    print("Loading DeepLoc IEDB predictions...")
    deeploc_iedb = pd.read_csv(
        deeploc_dir / "deeploc_IEDB_prediction.csv",
        dtype=str
    )

    print("Loading DeepLoc pathogen predictions...")
    deeploc_pathogen_files = sorted(
        deeploc_dir.glob("deeplocpro_prediction_*.csv")
    )

    deeploc_pathogen = pd.concat(
        [pd.read_csv(file, dtype=str) for file in deeploc_pathogen_files],
        ignore_index=True
    )

    return perfect_match, netsurfp_iedb, netsurfp_pathogen, deeploc_iedb, deeploc_pathogen


def clean_prediction_ids(
    perfect_match: pd.DataFrame,
    netsurfp_iedb: pd.DataFrame,
    netsurfp_pathogen: pd.DataFrame,
    deeploc_iedb: pd.DataFrame,
    deeploc_pathogen: pd.DataFrame,
):
    perfect_match = perfect_match.copy()
    netsurfp_iedb = netsurfp_iedb.copy()
    netsurfp_pathogen = netsurfp_pathogen.copy()
    deeploc_iedb = deeploc_iedb.copy()
    deeploc_pathogen = deeploc_pathogen.copy()

    for df in [deeploc_pathogen, deeploc_iedb, netsurfp_pathogen, netsurfp_iedb]:
        df.columns = df.columns.str.strip()
        df.drop(
            columns=df.columns[df.columns.str.contains("^Unnamed")],
            inplace=True
        )

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

    return perfect_match, netsurfp_iedb, netsurfp_pathogen, deeploc_iedb, deeploc_pathogen


def check_prediction_coverage(
    perfect_match: pd.DataFrame,
    netsurfp_iedb: pd.DataFrame,
    netsurfp_pathogen: pd.DataFrame,
    deeploc_iedb: pd.DataFrame,
    deeploc_pathogen: pd.DataFrame,
) -> None:
    matched_iedb_ids = set(
        perfect_match.loc[
            perfect_match["Matched"],
            "IEDB_Protein_ID_clean"
        ].dropna()
    )

    matched_pathogen_ids = set(
        perfect_match.loc[
            perfect_match["Matched"],
            "Pathogen_Protein_ID_clean"
        ].dropna()
    )

    missing_netsurfp_iedb = matched_iedb_ids - set(netsurfp_iedb["Protein_ID_clean"].dropna())
    missing_netsurfp_pathogen = matched_pathogen_ids - set(netsurfp_pathogen["Protein_ID_clean"].dropna())

    missing_deeploc_iedb = matched_iedb_ids - set(deeploc_iedb["Protein_ID_clean"].dropna())
    missing_deeploc_pathogen = matched_pathogen_ids - set(deeploc_pathogen["Protein_ID_clean"].dropna())

    print("Expected matched IEDB proteins:", len(matched_iedb_ids))
    print("Expected matched pathogen proteins:", len(matched_pathogen_ids))

    print("\nMissing from NetSurfP IEDB:", len(missing_netsurfp_iedb))
    print("Missing from NetSurfP pathogen:", len(missing_netsurfp_pathogen))

    print("\nMissing from DeepLoc IEDB:", len(missing_deeploc_iedb))
    print("Missing from DeepLoc pathogen:", len(missing_deeploc_pathogen))

    missing_netsurfp_pathogen_df = pd.DataFrame({
        "Pathogen_Protein_ID": sorted(missing_netsurfp_pathogen)
    })

    save_csv(
        missing_netsurfp_pathogen_df,
        DATA_DIR / "proccesed/missing_netsurfp_pathogen_ids.csv"
    )


def add_deeploc_metadata(
    perfect_match: pd.DataFrame,
    deeploc_iedb: pd.DataFrame,
    deeploc_pathogen: pd.DataFrame,
) -> pd.DataFrame:
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

    return perfect_match_meta


def prepare_numeric_columns(
    perfect_match_meta: pd.DataFrame,
    netsurfp_iedb: pd.DataFrame,
    netsurfp_pathogen: pd.DataFrame,
):
    for df in [netsurfp_iedb, netsurfp_pathogen]:
        df["n"] = pd.to_numeric(df["n"], errors="coerce")
        df["rsa"] = pd.to_numeric(df["rsa"], errors="coerce")

    for col in ["Pathogen_Start", "Pathogen_End", "Epitope_Start", "Epitope_End"]:
        perfect_match_meta[col] = pd.to_numeric(
            perfect_match_meta[col],
            errors="coerce"
        )

    return perfect_match_meta, netsurfp_iedb, netsurfp_pathogen


def add_rsa_metadata(
    perfect_match_meta: pd.DataFrame,
    netsurfp_iedb: pd.DataFrame,
    netsurfp_pathogen: pd.DataFrame,
) -> pd.DataFrame:
    perfect_match_meta = perfect_match_meta.reset_index(drop=True)
    perfect_match_meta["match_row_id"] = perfect_match_meta.index

    pathogen_regions = perfect_match_meta[
        [
            "match_row_id",
            "Pathogen_Protein_ID_clean",
            "Pathogen_Start",
            "Pathogen_End",
        ]
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

    iedb_regions = perfect_match_meta[
        [
            "match_row_id",
            "IEDB_Protein_ID_clean",
            "Epitope_Start",
            "Epitope_End",
        ]
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

    perfect_match_meta = (
        perfect_match_meta
        .merge(pathogen_rsa_summary, how="left", on="match_row_id")
        .merge(iedb_rsa_summary, how="left", on="match_row_id")
        .drop(columns=["match_row_id"])
    )

    return perfect_match_meta


def print_final_diagnostics(perfect_match_meta: pd.DataFrame) -> None:
    print("Rows:", len(perfect_match_meta))

    print("\nDeepLoc merge coverage:")
    print(
        "Pathogen location missing:",
        perfect_match_meta["pathogen_prot_location"].isna().sum()
    )
    print(
        "IEDB location missing:",
        perfect_match_meta["IEDB_prot_location"].isna().sum()
    )

    print("\nNetSurfP RSA coverage:")
    print(
        "Pathogen RSA missing:",
        perfect_match_meta["rsa_pathogen_peptide_mean"].isna().sum()
    )
    print(
        "IEDB RSA missing:",
        perfect_match_meta["rsa_epitope_mean"].isna().sum()
    )

    print(perfect_match_meta["Matched"].value_counts(dropna=False))


def main() -> None:
    (
        perfect_match,
        netsurfp_iedb,
        netsurfp_pathogen,
        deeploc_iedb,
        deeploc_pathogen,
    ) = load_prediction_inputs()

    (
        perfect_match,
        netsurfp_iedb,
        netsurfp_pathogen,
        deeploc_iedb,
        deeploc_pathogen,
    ) = clean_prediction_ids(
        perfect_match,
        netsurfp_iedb,
        netsurfp_pathogen,
        deeploc_iedb,
        deeploc_pathogen,
    )

    check_prediction_coverage(
        perfect_match,
        netsurfp_iedb,
        netsurfp_pathogen,
        deeploc_iedb,
        deeploc_pathogen,
    )

    perfect_match_meta = add_deeploc_metadata(
        perfect_match,
        deeploc_iedb,
        deeploc_pathogen,
    )

    perfect_match_meta, netsurfp_iedb, netsurfp_pathogen = prepare_numeric_columns(
        perfect_match_meta,
        netsurfp_iedb,
        netsurfp_pathogen,
    )

    perfect_match_meta = add_rsa_metadata(
        perfect_match_meta,
        netsurfp_iedb,
        netsurfp_pathogen,
    )

    print_final_diagnostics(perfect_match_meta)

    output_path = DATA_DIR / "proccesed/perfect_match_DN_add.csv"
    save_csv(perfect_match_meta, output_path)

    print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()