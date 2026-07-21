# Imports
import pandas as pd
from crossreactivity.io import save_csv, DATA_DIR
from crossreactivity.config import config


def standardize_iedb_columns(df: pd.DataFrame) -> pd.DataFrame:
    rename_dict = {
        "Assay ID - IEDB IRI": "Assay_ID",
        "Epitope - Name": "Sequence",
        "Epitope - Molecule Parent": "Protein_source",
        "Epitope - Molecule Parent IRI": "Protein_ID",
        "1st in vivo Process - Disease": "Disease",
        "1st in vivo Process - Process Type": "Disease_info",
        "1st in vivo Process - Disease Stage": "Disease_stage",
        "Epitope - Starting Position": "epitope_start_pos",
        "Epitope - Ending Position": "epitope_end_pos",
        "Epitope - Modified residues": "Modified_residues",
        "Assay - Response measured": "Response_measured",
        "Effector Cell - Name": "Effector_cell",
        "MHC Restriction - Name": "MHC_restriction"
    }

    cols_to_keep = list(rename_dict.values())

    return (
        df
        .rename(columns=rename_dict)
        .filter(cols_to_keep)
    )


def clean_iedb_data(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    df["Disease"] = (
        df["Disease"]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    df["Sequence"] = (
        df["Sequence"]
        .fillna("")
        .astype(str)
        .str.replace(r"\+.*", "", regex=True)
        .str.strip()
    )

    df["Protein_ID"] = (
        df["Protein_ID"]
        .fillna("")
        .astype(str)
        .str.extract(r"([^/]+$)", expand=False)
    )

    return df


def remove_nested_epitopes(df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove nested epitopes, i.e. sequences that are fully contained
    within another sequence from the same protein.
    """

    grouped_results = []

    for protein_id, group in df.groupby("Protein_ID", sort=False):
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

    df = pd.concat(grouped_results, ignore_index=True)

    return df[~df["Is_nested"]].drop(columns=["Is_nested"])


def main() -> None:
    iedb_df = pd.read_csv(DATA_DIR / "raw" / "tcell_table_export_1781203145.csv")

    print(f"IEDB dataset size: {len(iedb_df)}")

    # Standardize and clean columns
    iedb_df = standardize_iedb_columns(iedb_df)
    iedb_df = clean_iedb_data(iedb_df)

    print(f"Dataset size before removing non-autoimmune: {len(iedb_df)}")
    print(f"Unique all assays: {iedb_df['Assay_ID'].nunique()}")

    # ------------------------------------------------------------
    # 1. Keep autoimmune rows only
    # ------------------------------------------------------------
    autoimmune = (
        iedb_df[iedb_df["Disease_info"] == "Occurrence of autoimmune disease"]
        .copy()
        .drop(columns=["Disease_info"])
    )

    print(f"Dataset size after removing non-autoimmune: {len(autoimmune)}")
    print(f"Unique autoimmune assays: {autoimmune['Assay_ID'].nunique()}")

    # ------------------------------------------------------------
    # 2. Drop duplicate Protein_ID + Sequence + Modified_residues
    # ------------------------------------------------------------
    print(f"Dataset size before duplicate filtering: {len(autoimmune)}")
    print(f"Unique assays before duplicate filtering: {autoimmune['Assay_ID'].nunique()}")

    unique_epitopes = autoimmune.drop_duplicates(
        subset=["Protein_ID", "Sequence", "Modified_residues"],
        keep="first"
    ).copy()

    print(f"Dataset size after duplicate filtering: {len(unique_epitopes)}")
    print(f"Unique assays after duplicate filtering: {unique_epitopes['Assay_ID'].nunique()}")
    print(
        "Unique Protein_ID + Sequence + Modified_residues groups:",
        unique_epitopes[["Protein_ID", "Sequence", "Modified_residues"]]
        .drop_duplicates()
        .shape[0]
    )

    # Get length limits from config
    min_len = config["filtering"]["min_epitope_length"]
    max_len = config["filtering"]["max_epitope_length"]

    # ------------------------------------------------------------
    # 3. Remove modified residues
    # ------------------------------------------------------------
    unmodified = unique_epitopes[
        unique_epitopes["Modified_residues"].isna()
    ].copy()

    print(f"Dataset size after filtering modified residues: {len(unmodified)}")
    print(f"Unique assays after filtering modified residues: {unmodified['Assay_ID'].nunique()}")

    # ------------------------------------------------------------
    # 4. Filter epitope length
    # ------------------------------------------------------------
    length_filtered = unmodified[
        unmodified["Sequence"].str.len().between(
            min_len,
            max_len,
            inclusive="both"
        )
    ].copy()

    print(f"Dataset size after filtering {min_len}-{max_len}: {len(length_filtered)}")
    print(f"Unique assays after length filtering: {length_filtered['Assay_ID'].nunique()}")

    # ------------------------------------------------------------
    # 5. Remove nested epitopes
    # ------------------------------------------------------------
    filtered_sequences = remove_nested_epitopes(length_filtered)

    print(f"Dataset size after removing nested epitopes: {len(filtered_sequences)}")
    print(f"Unique retained assays after nested filtering: {filtered_sequences['Assay_ID'].nunique()}")
    print(f"Unique source proteins after final filtering: {filtered_sequences['Protein_ID'].nunique()}")
    print(
        "Unique final protein-sequence pairs:",
        filtered_sequences[["Protein_ID", "Sequence"]].drop_duplicates().shape[0]
    )

    # ------------------------------------------------------------
    # 6. Save
    # ------------------------------------------------------------
    output_path = DATA_DIR / "intermediate/wrangled_IEDB.csv"
    save_csv(filtered_sequences, output_path)

    print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()