import pandas as pd

from crossreactivity.io import DATA_DIR, load_csv, save_csv
from crossreactivity.reference_data import GRAM_STATUS


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


def concat_csv_files(files, description: str) -> pd.DataFrame:
    if not files:
        raise FileNotFoundError(f"No files found for: {description}")

    return pd.concat(
        [pd.read_csv(file, dtype=str) for file in files],
        ignore_index=True
    )


def concat_netsurfp_files(files, description: str) -> pd.DataFrame:
    if not files:
        raise FileNotFoundError(f"No files found for: {description}")

    return pd.concat(
        [read_netsurfp_file(file) for file in files],
        ignore_index=True
    )


def load_prediction_inputs():
    netsurfp_dir = DATA_DIR / "raw" / "netsurfp"
    deeploc_dir = DATA_DIR / "raw" / "deeploc"

    print("Loading region-labeled long match results...")
    match_df = load_csv(DATA_DIR / "proccesed/iedb_match_regions_long.csv")

    print("Loading pathogen protein metadata...")
    protein_sequences = load_csv(
        DATA_DIR / "intermediate/protein_sequences.csv"
    )

    print("Loading pathogen DeepLoc prediction manifest...")
    deeploc_manifest = load_csv(
        DATA_DIR / "intermediate/pathogen_deeploc_prediction_manifest.csv"
    )

    print("Loading IEDB metadata...")
    iedb_data = load_csv(DATA_DIR / "intermediate/wrangled_IEDB.csv")

    print("Loading NetSurfP IEDB predictions...")
    netsurfp_iedb_files = sorted(
        netsurfp_dir.glob("netsurfp_IEDB_prediction_*.csv")
    )

    netsurfp_iedb = concat_netsurfp_files(
        netsurfp_iedb_files,
        "NetSurfP IEDB predictions"
    )

    print("Loading NetSurfP pathogen predictions...")
    netsurfp_pathogen_files = sorted(
        netsurfp_dir.glob("netsurfp_pathogen_prediction_*.csv")
    )

    netsurfp_pathogen = concat_netsurfp_files(
        netsurfp_pathogen_files,
        "NetSurfP pathogen predictions"
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

    deeploc_pathogen = concat_csv_files(
        deeploc_pathogen_files,
        "DeepLoc pathogen predictions"
    )

    return (
        match_df,
        protein_sequences,
        deeploc_manifest,
        iedb_data,
        netsurfp_iedb,
        netsurfp_pathogen,
        deeploc_iedb,
        deeploc_pathogen,
    )


def clean_prediction_ids(
    match_df: pd.DataFrame,
    protein_sequences: pd.DataFrame,
    deeploc_manifest: pd.DataFrame,
    iedb_data: pd.DataFrame,
    netsurfp_iedb: pd.DataFrame,
    netsurfp_pathogen: pd.DataFrame,
    deeploc_iedb: pd.DataFrame,
    deeploc_pathogen: pd.DataFrame,
):
    match_df = match_df.copy()
    protein_sequences = protein_sequences.copy()
    deeploc_manifest = deeploc_manifest.copy()
    iedb_data = iedb_data.copy()
    netsurfp_iedb = netsurfp_iedb.copy()
    netsurfp_pathogen = netsurfp_pathogen.copy()
    deeploc_iedb = deeploc_iedb.copy()
    deeploc_pathogen = deeploc_pathogen.copy()

    for df in [
        match_df,
        protein_sequences,
        deeploc_manifest,
        iedb_data,
        netsurfp_iedb,
        netsurfp_pathogen,
        deeploc_iedb,
        deeploc_pathogen,
    ]:
        df.columns = df.columns.str.strip()
        df.drop(
            columns=df.columns[df.columns.str.contains("^Unnamed")],
            inplace=True
        )

    match_df["Pathogen_Protein_ID_clean"] = (
        match_df["Pathogen_Protein_ID"]
        .astype("string")
        .str.strip()
        .str.upper()
    )

    match_df["IEDB_Protein_ID_clean"] = (
        match_df["IEDB_Protein_ID"]
        .astype("string")
        .str.strip()
        .str.upper()
    )

    protein_sequences["Protein_ID_clean"] = (
        protein_sequences["Protein_ID"]
        .astype("string")
        .str.strip()
        .str.upper()
    )

    if protein_sequences["Protein_ID_clean"].isna().any():
        raise ValueError("protein_sequences.csv contains missing Protein_IDs.")

    if protein_sequences["Protein_ID_clean"].duplicated().any():
        raise ValueError(
            "protein_sequences.csv must contain one row per Protein_ID."
        )

    required_manifest_cols = {"Prediction_ID", "Protein_ID", "Gram_status"}
    missing_manifest_cols = required_manifest_cols - set(deeploc_manifest.columns)
    if missing_manifest_cols:
        raise ValueError(
            "DeepLoc manifest is missing columns: "
            + ", ".join(sorted(missing_manifest_cols))
        )

    deeploc_manifest["Prediction_ID_clean"] = (
        deeploc_manifest["Prediction_ID"]
        .astype("string")
        .str.strip()
        .str.upper()
    )
    deeploc_manifest["Protein_ID_clean"] = (
        deeploc_manifest["Protein_ID"]
        .astype("string")
        .str.strip()
        .str.upper()
    )
    deeploc_manifest["Gram_status"] = (
        deeploc_manifest["Gram_status"]
        .astype("string")
        .str.strip()
        .str.lower()
    )

    manifest_keys = ["Protein_ID_clean", "Gram_status"]
    if deeploc_manifest[manifest_keys].isna().any(axis=None):
        raise ValueError("DeepLoc manifest contains missing protein or Gram keys.")
    if deeploc_manifest.duplicated(manifest_keys).any():
        raise ValueError(
            "DeepLoc manifest has duplicate Protein_ID + Gram_status rows."
        )
    if deeploc_manifest["Prediction_ID_clean"].duplicated().any():
        raise ValueError("DeepLoc manifest contains duplicate Prediction_IDs.")

    if "Pathogen_Organism" not in match_df.columns:
        raise ValueError("Match table is missing Pathogen_Organism.")

    match_df["Gram_status"] = (
        match_df["Pathogen_Organism"]
        .astype("string")
        .str.strip()
        .map(GRAM_STATUS)
    )

    unknown_organisms = (
        match_df.loc[match_df["Gram_status"].isna(), "Pathogen_Organism"]
        .drop_duplicates()
        .sort_values()
    )
    if not unknown_organisms.empty:
        raise ValueError(
            "Missing Gram-status mapping for: "
            + ", ".join(unknown_organisms.astype(str))
        )

    iedb_data["Protein_ID_clean"] = (
        iedb_data["Protein_ID"]
        .astype("string")
        .str.strip()
        .str.upper()
    )

    iedb_data["Assay_ID"] = (
        iedb_data["Assay_ID"]
        .astype("string")
        .str.strip()
    )

    match_df["Assay_ID"] = (
        match_df["Assay_ID"]
        .astype("string")
        .str.strip()
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

    deeploc_pathogen["Prediction_ID_clean"] = (
        deeploc_pathogen["ACC"]
        .astype("string")
        .str.strip()
        .str.upper()
    )

    return (
        match_df,
        protein_sequences,
        deeploc_manifest,
        iedb_data,
        netsurfp_iedb,
        netsurfp_pathogen,
        deeploc_iedb,
        deeploc_pathogen,
    )


def add_input_metadata(
    match_df: pd.DataFrame,
    protein_sequences: pd.DataFrame,
    iedb_data: pd.DataFrame,
) -> pd.DataFrame:
    """
    Add metadata that is not already present in the match table.
    """
    match_df = match_df.copy()
    protein_sequences = protein_sequences.copy()
    iedb_data = iedb_data.copy()

    # ------------------------------------------------------------
    # Pathogen protein length
    # ------------------------------------------------------------
    pathogen_meta = protein_sequences[
        [
            "Protein_ID_clean",
            "Sequence",
        ]
    ].copy()

    pathogen_meta["Pathogen_Protein_Length"] = (
        pathogen_meta["Sequence"]
        .astype("string")
        .str.len()
    )

    pathogen_meta = (
        pathogen_meta
        .drop(columns="Sequence")
    )

    before_rows = len(match_df)

    match_df = (
        match_df
        .merge(
            pathogen_meta,
            how="left",
            left_on="Pathogen_Protein_ID_clean",
            right_on="Protein_ID_clean",
            validate="many_to_one",
        )
        .drop(columns="Protein_ID_clean")
    )

    if len(match_df) != before_rows:
        raise RuntimeError(
            "Pathogen metadata merge changed the number of match rows."
        )

    if match_df["Pathogen_Protein_Length"].isna().any():
        missing = (
            match_df.loc[
                match_df["Pathogen_Protein_Length"].isna(),
                "Pathogen_Protein_ID_clean",
            ]
            .dropna()
            .unique()[:10]
        )
        raise ValueError(
            "Matched pathogen proteins missing from protein_sequences.csv: "
            + ", ".join(missing)
        )

    # ------------------------------------------------------------
    # IEDB MHC restriction
    # ------------------------------------------------------------
    iedb_meta = (
        iedb_data[
            [
                "Assay_ID",
                "MHC_restriction",
            ]
        ]
        .drop_duplicates(subset="Assay_ID")
    )

    before_rows = len(match_df)

    match_df = match_df.merge(
        iedb_meta,
        how="left",
        on="Assay_ID",
        validate="many_to_one",
    )

    if len(match_df) != before_rows:
        raise RuntimeError(
            "IEDB metadata merge changed the number of match rows."
        )

    return match_df


def check_prediction_coverage(
    match_df: pd.DataFrame,
    netsurfp_iedb: pd.DataFrame,
    netsurfp_pathogen: pd.DataFrame,
    deeploc_iedb: pd.DataFrame,
    deeploc_pathogen: pd.DataFrame,
) -> None:
    matched_iedb_ids = set(
        match_df["IEDB_Protein_ID_clean"]
        .dropna()
        .tolist()
    )

    matched_pathogen_ids = set(
        match_df["Pathogen_Protein_ID_clean"]
        .dropna()
        .tolist()
    )

    expected_pathogen_prediction_ids = set(
        match_df["Prediction_ID_clean"]
        .dropna()
        .tolist()
    )

    missing_netsurfp_iedb = (
        matched_iedb_ids -
        set(netsurfp_iedb["Protein_ID_clean"].dropna())
    )

    missing_netsurfp_pathogen = (
        matched_pathogen_ids -
        set(netsurfp_pathogen["Protein_ID_clean"].dropna())
    )

    missing_deeploc_iedb = (
        matched_iedb_ids -
        set(deeploc_iedb["Protein_ID_clean"].dropna())
    )

    missing_deeploc_pathogen = (
        expected_pathogen_prediction_ids -
        set(deeploc_pathogen["Prediction_ID_clean"].dropna())
    )

    print("Expected matched IEDB proteins:", len(matched_iedb_ids))
    print("Expected matched pathogen proteins:", len(matched_pathogen_ids))
    print(
        "Expected Gram-specific pathogen DeepLoc predictions:",
        len(expected_pathogen_prediction_ids),
    )

    print("\nMissing from NetSurfP IEDB:", len(missing_netsurfp_iedb))
    print("Missing from NetSurfP pathogen:", len(missing_netsurfp_pathogen))

    print("\nMissing from DeepLoc IEDB:", len(missing_deeploc_iedb))
    print("Missing from DeepLoc pathogen:", len(missing_deeploc_pathogen))

    save_csv(
        pd.DataFrame({
            "IEDB_Protein_ID": sorted(missing_netsurfp_iedb)
        }),
        DATA_DIR / "proccesed/missing_netsurfp_iedb_ids.csv"
    )

    save_csv(
        pd.DataFrame({
            "Pathogen_Protein_ID": sorted(missing_netsurfp_pathogen)
        }),
        DATA_DIR / "proccesed/missing_netsurfp_pathogen_ids.csv"
    )

    save_csv(
        pd.DataFrame({
            "IEDB_Protein_ID": sorted(missing_deeploc_iedb)
        }),
        DATA_DIR / "proccesed/missing_deeploc_iedb_ids.csv"
    )

    missing_deeploc_pathogen_df = (
        match_df.loc[
            match_df["Prediction_ID_clean"].isin(missing_deeploc_pathogen),
            [
                "Prediction_ID_clean",
                "Pathogen_Protein_ID_clean",
                "Gram_status",
            ],
        ]
        .drop_duplicates()
        .rename(columns={
            "Prediction_ID_clean": "Prediction_ID",
            "Pathogen_Protein_ID_clean": "Protein_ID",
        })
        .sort_values("Prediction_ID")
    )
    save_csv(
        missing_deeploc_pathogen_df,
        DATA_DIR / "proccesed/missing_deeploc_pathogen_ids.csv"
    )


def add_pathogen_prediction_ids(
    match_df: pd.DataFrame,
    deeploc_manifest: pd.DataFrame,
) -> pd.DataFrame:
    before_rows = len(match_df)
    manifest = deeploc_manifest[
        ["Protein_ID_clean", "Gram_status", "Prediction_ID_clean"]
    ]

    match_df = match_df.merge(
        manifest,
        how="left",
        left_on=["Pathogen_Protein_ID_clean", "Gram_status"],
        right_on=["Protein_ID_clean", "Gram_status"],
        validate="many_to_one",
    ).drop(columns="Protein_ID_clean")

    if len(match_df) != before_rows:
        raise RuntimeError("DeepLoc manifest merge changed match row count.")

    missing = match_df["Prediction_ID_clean"].isna()
    if missing.any():
        examples = (
            match_df.loc[
                missing,
                ["Pathogen_Protein_ID_clean", "Gram_status"],
            ]
            .drop_duplicates()
            .head(10)
            .apply(lambda row: f"{row.iloc[0]} ({row.iloc[1]})", axis=1)
        )
        raise ValueError(
            "Matched protein/Gram pairs missing from the DeepLoc manifest: "
            + ", ".join(examples)
        )

    return match_df


def unique_prediction_metadata(
    df: pd.DataFrame,
    key: str,
    value: str,
    description: str,
) -> pd.DataFrame:
    metadata = df[[key, value]].dropna(subset=[key]).drop_duplicates()
    conflicts = metadata.groupby(key, dropna=False)[value].nunique(dropna=False)
    conflicts = conflicts[conflicts > 1]

    if not conflicts.empty:
        examples = ", ".join(conflicts.index.astype(str)[:10])
        raise ValueError(
            f"Conflicting {description} predictions for: {examples}"
        )

    return metadata.drop_duplicates(key)


def add_deeploc_metadata(
    match_df: pd.DataFrame,
    deeploc_iedb: pd.DataFrame,
    deeploc_pathogen: pd.DataFrame,
) -> pd.DataFrame:
    iedb_location_meta = unique_prediction_metadata(
        deeploc_iedb,
        key="Protein_ID_clean",
        value="Localizations",
        description="IEDB DeepLoc",
    ).rename(
        columns={"Localizations": "IEDB_prot_location"}
    )

    pathogen_location_meta = unique_prediction_metadata(
        deeploc_pathogen,
        key="Prediction_ID_clean",
        value="Localization",
        description="pathogen DeepLoc",
    ).rename(
        columns={"Localization": "pathogen_prot_location"}
    )

    before_rows = len(match_df)
    match_df_meta = match_df.merge(
        pathogen_location_meta,
        how="left",
        on="Prediction_ID_clean",
        validate="many_to_one",
    )

    match_df_meta = (
        match_df_meta.merge(
            iedb_location_meta,
            how="left",
            left_on="IEDB_Protein_ID_clean",
            right_on="Protein_ID_clean",
            validate="many_to_one",
        )
        .drop(columns=["Protein_ID_clean", "Gram_status", "Prediction_ID_clean"])
    )

    if len(match_df_meta) != before_rows:
        raise RuntimeError("DeepLoc merges changed the number of match rows.")

    return match_df_meta


def prepare_numeric_columns(
    match_df_meta: pd.DataFrame,
    netsurfp_iedb: pd.DataFrame,
    netsurfp_pathogen: pd.DataFrame,
):
    for df in [netsurfp_iedb, netsurfp_pathogen]:
        df["n"] = pd.to_numeric(df["n"], errors="coerce")
        df["rsa"] = pd.to_numeric(df["rsa"], errors="coerce")

    coordinate_cols = [
        "Pathogen_Start",
        "Pathogen_End",
        "Epitope_Start",
        "Epitope_End",
        "IEDB_Match_Start",
        "IEDB_Match_End",
        "IEDB_Region_Start",
        "IEDB_Region_End",
    ]

    for col in coordinate_cols:
        if col in match_df_meta.columns:
            match_df_meta[col] = pd.to_numeric(
                match_df_meta[col],
                errors="coerce"
            )

    if "Pathogen_Protein_Length" in match_df_meta.columns:
        match_df_meta["Pathogen_Protein_Length"] = pd.to_numeric(
            match_df_meta["Pathogen_Protein_Length"],
            errors="coerce"
        )

    return match_df_meta, netsurfp_iedb, netsurfp_pathogen


def build_rsa_lookup(netsurfp_df: pd.DataFrame) -> dict:
    """
    Build lookup:
        Protein_ID_clean -> residue-level RSA table

    This avoids huge row expansion from merging every match row
    with every residue row for the same protein.
    """
    netsurfp_df = netsurfp_df.dropna(
        subset=["Protein_ID_clean", "n", "rsa"]
    ).copy()

    rsa_lookup = {}

    for protein_id, group in netsurfp_df.groupby("Protein_ID_clean"):
        rsa_lookup[protein_id] = (
            group[["n", "rsa"]]
            .sort_values("n")
            .reset_index(drop=True)
        )

    return rsa_lookup


def mean_rsa_for_region(
    rsa_lookup: dict,
    protein_id,
    start,
    end,
):
    if pd.isna(protein_id) or pd.isna(start) or pd.isna(end):
        return None

    protein_id = str(protein_id).strip().upper()

    if protein_id not in rsa_lookup:
        return None

    try:
        start = int(start)
        end = int(end)
    except ValueError:
        return None

    if start > end:
        start, end = end, start

    rsa_df = rsa_lookup[protein_id]

    region_rsa = rsa_df.loc[
        rsa_df["n"].between(start, end),
        "rsa"
    ]

    if region_rsa.empty:
        return None

    return region_rsa.mean()


def add_rsa_metadata(
    match_df_meta: pd.DataFrame,
    netsurfp_iedb: pd.DataFrame,
    netsurfp_pathogen: pd.DataFrame,
) -> pd.DataFrame:
    match_df_meta = match_df_meta.copy()

    print("Building NetSurfP RSA lookup tables...")

    pathogen_rsa_lookup = build_rsa_lookup(netsurfp_pathogen)
    iedb_rsa_lookup = build_rsa_lookup(netsurfp_iedb)

    print("Adding pathogen matched-peptide RSA means...")

    match_df_meta["rsa_pathogen_peptide_mean"] = [
        mean_rsa_for_region(
            rsa_lookup=pathogen_rsa_lookup,
            protein_id=protein_id,
            start=start,
            end=end,
        )
        for protein_id, start, end in zip(
            match_df_meta["Pathogen_Protein_ID_clean"],
            match_df_meta["Pathogen_Start"],
            match_df_meta["Pathogen_End"],
        )
    ]

    print("Adding IEDB matched-region RSA means...")

    match_df_meta["rsa_iedb_matched_region_mean"] = [
        mean_rsa_for_region(
            rsa_lookup=iedb_rsa_lookup,
            protein_id=protein_id,
            start=start,
            end=end,
        )
        for protein_id, start, end in zip(
            match_df_meta["IEDB_Protein_ID_clean"],
            match_df_meta["IEDB_Match_Start"],
            match_df_meta["IEDB_Match_End"],
        )
    ]

    print("Adding full IEDB epitope RSA means...")

    match_df_meta["rsa_epitope_mean"] = [
        mean_rsa_for_region(
            rsa_lookup=iedb_rsa_lookup,
            protein_id=protein_id,
            start=start,
            end=end,
        )
        for protein_id, start, end in zip(
            match_df_meta["IEDB_Protein_ID_clean"],
            match_df_meta["Epitope_Start"],
            match_df_meta["Epitope_End"],
        )
    ]

    return match_df_meta


def print_final_diagnostics(match_df_meta: pd.DataFrame) -> None:
    print("Rows:", len(match_df_meta))

    print("\nUnique IDs:")
    print(
        "Unique pathogen proteins:",
        match_df_meta["Pathogen_Protein_ID_clean"].nunique()
    )
    print(
        "Unique IEDB proteins:",
        match_df_meta["IEDB_Protein_ID_clean"].nunique()
    )

    if "Assay_ID" in match_df_meta.columns:
        print("Unique assays:", match_df_meta["Assay_ID"].nunique())

    if "IEDB_Region_ID" in match_df_meta.columns:
        print("Unique IEDB regions:", match_df_meta["IEDB_Region_ID"].nunique())

    print("\nInput metadata merge coverage:")

    if "Pathogen_Organism" in match_df_meta.columns:
        print(
            "Pathogen organism missing:",
            match_df_meta["Pathogen_Organism"].isna().sum()
        )

    if "Pathogen_Protein_Length" in match_df_meta.columns:
        print(
            "Pathogen protein length missing:",
            match_df_meta["Pathogen_Protein_Length"].isna().sum()
        )

    if "MHC_restriction" in match_df_meta.columns:
        print(
            "IEDB MHC restriction missing:",
            match_df_meta["MHC_restriction"].isna().sum()
        )

    print("\nDeepLoc merge coverage:")
    print(
        "Pathogen location missing:",
        match_df_meta["pathogen_prot_location"].isna().sum()
    )
    print(
        "IEDB location missing:",
        match_df_meta["IEDB_prot_location"].isna().sum()
    )

    print("\nNetSurfP RSA coverage:")
    print(
        "Pathogen matched peptide RSA missing:",
        match_df_meta["rsa_pathogen_peptide_mean"].isna().sum()
    )
    print(
        "IEDB matched region RSA missing:",
        match_df_meta["rsa_iedb_matched_region_mean"].isna().sum()
    )
    print(
        "IEDB full epitope RSA missing:",
        match_df_meta["rsa_epitope_mean"].isna().sum()
    )


def main() -> None:
    (
        match_df,
        protein_sequences,
        deeploc_manifest,
        iedb_data,
        netsurfp_iedb,
        netsurfp_pathogen,
        deeploc_iedb,
        deeploc_pathogen,
    ) = load_prediction_inputs()

    (
        match_df,
        protein_sequences,
        deeploc_manifest,
        iedb_data,
        netsurfp_iedb,
        netsurfp_pathogen,
        deeploc_iedb,
        deeploc_pathogen,
    ) = clean_prediction_ids(
        match_df,
        protein_sequences,
        deeploc_manifest,
        iedb_data,
        netsurfp_iedb,
        netsurfp_pathogen,
        deeploc_iedb,
        deeploc_pathogen,
    )

    match_df = add_input_metadata(
        match_df=match_df,
        protein_sequences=protein_sequences,
        iedb_data=iedb_data,
    )

    match_df = add_pathogen_prediction_ids(
        match_df=match_df,
        deeploc_manifest=deeploc_manifest,
    )

    check_prediction_coverage(
        match_df,
        netsurfp_iedb,
        netsurfp_pathogen,
        deeploc_iedb,
        deeploc_pathogen,
    )

    match_df_meta = add_deeploc_metadata(
        match_df,
        deeploc_iedb,
        deeploc_pathogen,
    )

    match_df_meta, netsurfp_iedb, netsurfp_pathogen = prepare_numeric_columns(
        match_df_meta,
        netsurfp_iedb,
        netsurfp_pathogen,
    )

    match_df_meta = add_rsa_metadata(
        match_df_meta,
        netsurfp_iedb,
        netsurfp_pathogen,
    )

    print_final_diagnostics(match_df_meta)

    output_path = (DATA_DIR / "proccesed/iedb_match_regions_long_metadata.csv")
    save_csv(match_df_meta, output_path)

    print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()