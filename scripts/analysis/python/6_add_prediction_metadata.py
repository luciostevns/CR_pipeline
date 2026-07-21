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
    protein_meta = load_csv(DATA_DIR / "intermediate/protein_metadata.csv")

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
        protein_meta,
        iedb_data,
        netsurfp_iedb,
        netsurfp_pathogen,
        deeploc_iedb,
        deeploc_pathogen,
    )


def clean_prediction_ids(
    match_df: pd.DataFrame,
    protein_meta: pd.DataFrame,
    iedb_data: pd.DataFrame,
    netsurfp_iedb: pd.DataFrame,
    netsurfp_pathogen: pd.DataFrame,
    deeploc_iedb: pd.DataFrame,
    deeploc_pathogen: pd.DataFrame,
):
    match_df = match_df.copy()
    protein_meta = protein_meta.copy()
    iedb_data = iedb_data.copy()
    netsurfp_iedb = netsurfp_iedb.copy()
    netsurfp_pathogen = netsurfp_pathogen.copy()
    deeploc_iedb = deeploc_iedb.copy()
    deeploc_pathogen = deeploc_pathogen.copy()

    for df in [
        match_df,
        protein_meta,
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

    if "Proteome_ID" in match_df.columns:
        match_df["Proteome_ID"] = (
            match_df["Proteome_ID"]
            .astype("string")
            .str.strip()
        )

    protein_meta["Protein_ID_clean"] = (
        protein_meta["Protein_ID"]
        .astype("string")
        .str.strip()
        .str.upper()
    )

    if "Proteome_ID" in protein_meta.columns:
        protein_meta["Proteome_ID"] = (
            protein_meta["Proteome_ID"]
            .astype("string")
            .str.strip()
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

    deeploc_pathogen["Protein_ID_clean"] = (
        deeploc_pathogen["ACC"]
        .astype(str)
        .str.strip()
        .str.upper()
    )

    return (
        match_df,
        protein_meta,
        iedb_data,
        netsurfp_iedb,
        netsurfp_pathogen,
        deeploc_iedb,
        deeploc_pathogen,
    )


def add_input_metadata(
    match_df: pd.DataFrame,
    protein_meta: pd.DataFrame,
    iedb_data: pd.DataFrame,
) -> pd.DataFrame:
    """
    Add metadata that is not already present in the match table.
    """
    match_df = match_df.copy()
    protein_meta = protein_meta.copy()
    iedb_data = iedb_data.copy()

    # ------------------------------------------------------------
    # Pathogen protein length
    # ------------------------------------------------------------
    pathogen_meta = protein_meta[
        [
            "Protein_ID_clean",
            "Proteome_ID",
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
        .drop_duplicates(
            subset=["Protein_ID_clean", "Proteome_ID"]
        )
    )

    before_rows = len(match_df)

    match_df = (
        match_df
        .merge(
            pathogen_meta,
            how="left",
            left_on=["Pathogen_Protein_ID_clean", "Proteome_ID"],
            right_on=["Protein_ID_clean", "Proteome_ID"],
            validate="many_to_one",
        )
        .drop(columns="Protein_ID_clean")
    )

    if len(match_df) != before_rows:
        raise RuntimeError(
            "Pathogen metadata merge changed the number of match rows."
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
        matched_pathogen_ids -
        set(deeploc_pathogen["Protein_ID_clean"].dropna())
    )

    print("Expected matched IEDB proteins:", len(matched_iedb_ids))
    print("Expected matched pathogen proteins:", len(matched_pathogen_ids))

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

    save_csv(
        pd.DataFrame({
            "Pathogen_Protein_ID": sorted(missing_deeploc_pathogen)
        }),
        DATA_DIR / "proccesed/missing_deeploc_pathogen_ids.csv"
    )


def add_deeploc_metadata(
    match_df: pd.DataFrame,
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

    match_df_meta = (
        match_df
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
            right_on="Protein_ID_clean")
        .drop(columns=["Protein_ID_clean"])
    )

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

    if "Pathogen_Metadata_Protein_Length" in match_df_meta.columns:
        match_df_meta["Pathogen_Metadata_Protein_Length"] = pd.to_numeric(
            match_df_meta["Pathogen_Metadata_Protein_Length"],
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

    if "Pathogen_Metadata_Genus_species" in match_df_meta.columns:
        print(
            "Pathogen genus/species missing:",
            match_df_meta["Pathogen_Metadata_Genus_species"].isna().sum()
        )

    if "Pathogen_Metadata_Protein_Length" in match_df_meta.columns:
        print(
            "Pathogen protein length missing:",
            match_df_meta["Pathogen_Metadata_Protein_Length"].isna().sum()
        )

    if "IEDB_Metadata_MHC_restriction" in match_df_meta.columns:
        print(
            "IEDB MHC restriction missing:",
            match_df_meta["IEDB_Metadata_MHC_restriction"].isna().sum()
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
        protein_meta,
        iedb_data,
        netsurfp_iedb,
        netsurfp_pathogen,
        deeploc_iedb,
        deeploc_pathogen,
    ) = load_prediction_inputs()

    (
        match_df,
        protein_meta,
        iedb_data,
        netsurfp_iedb,
        netsurfp_pathogen,
        deeploc_iedb,
        deeploc_pathogen,
    ) = clean_prediction_ids(
        match_df,
        protein_meta,
        iedb_data,
        netsurfp_iedb,
        netsurfp_pathogen,
        deeploc_iedb,
        deeploc_pathogen,
    )

    match_df = add_input_metadata(
        match_df=match_df,
        protein_meta=protein_meta,
        iedb_data=iedb_data,
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