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

    iedb_df = standardize_iedb_columns(iedb_df)

    iedb_df = clean_iedb_data(iedb_df)

    autoimmune = iedb_df[iedb_df["Disease_info"] == "Occurrence of autoimmune disease"].copy().drop(columns=["Disease_info"])

    print(f"Dataset size before removing non-autoimmune: {len(iedb_df)}")
    print(f"Unique all epitopes: {iedb_df['Assay_ID'].nunique()}")
    print(f"Dataset size after removing non-autoimmune: {len(autoimmune)}")
    print(f"Unique autoimmune epitopes: {autoimmune['Assay_ID'].nunique()}")

    autoimmune_unique = autoimmune.drop_duplicates(
        subset=["Protein_ID", "Sequence"]
    )

    print(f"Dataset size before removing duplicate seqs: {len(autoimmune)}")
    print(f"Dataset size after removing duplicate seqs: {len(autoimmune_unique)}")
    print(f"Unique autoimmune epitopes: {autoimmune_unique['Assay_ID'].nunique()}")

    min_len = config["filtering"]["min_epitope_length"]
    max_len = config["filtering"]["max_epitope_length"]

    mod_res_df = autoimmune_unique[autoimmune_unique["Modified_residues"].isna()].copy()

    print(f"Dataset size after filtering modified residues: {len(mod_res_df)}")

    filtered_data = mod_res_df[
        mod_res_df["Sequence"].str.len().between(min_len, max_len, inclusive="both")
    ].copy()

    print(f"Dataset size after filtering 12-25: {len(filtered_data)}")

    filtered_sequences = remove_nested_epitopes(filtered_data)

    print(f"Dataset size after removing nested epitopes: {len(filtered_sequences)}")

    output_path = DATA_DIR / "intermediate/wrangled_IEDB.csv"
    save_csv(filtered_sequences, output_path)

    print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()