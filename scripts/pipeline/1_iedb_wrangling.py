# Imports
import pandas as pd
from crossreactivity.io import load_excel, save_csv, DATA_DIR
from crossreactivity.config import config


def standardize_iedb_columns(df: pd.DataFrame) -> pd.DataFrame:
    rename_dict = {
        "Assay ID - IEDB IRI": "Assay_ID",
        "Epitope - Name": "Sequence",
        "Epitope - Molecule Parent": "Protein_source",
        "Epitope - Molecule Parent IRI": "Protein_ID",
        "1st in vivo Process - Disease": "Disease",
        "Epitope - Starting Position": "epitope_start_pos",
        "Epitope - Ending Position": "epitope_end_pos",
        "Epitope - Modified residues": "Modified_residues",
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
    autoimmune_ref = load_excel(DATA_DIR / "raw/IEDB_autoimmune_epitope_assays_raw.xlsx")
    general_iedb_data = load_excel(DATA_DIR / "raw/IEDB_epitop_assays_raw.xlsx")

    print(f"Autoimmune reference dataset size: {len(autoimmune_ref)}")
    print(f"General IEDB dataset size: {len(general_iedb_data)}")

    autoimmune_ref = standardize_iedb_columns(autoimmune_ref)
    general_iedb_data = standardize_iedb_columns(general_iedb_data)

    autoimmune_ref = clean_iedb_data(autoimmune_ref)
    general_iedb_data = clean_iedb_data(general_iedb_data)

    # Supplement autoimmune export with type 1 diabetes rows from the general IEDB export
    extra_data = general_iedb_data[
        general_iedb_data["Disease"].eq("type 1 diabetes mellitus")
    ].copy()

    print(f"Extra type 1 diabetes rows added: {len(extra_data)}")
    print(extra_data["Disease"].value_counts())

    combined_data = pd.concat([autoimmune_ref, extra_data], ignore_index=True)
    combined_data = combined_data.drop_duplicates(subset=["Assay_ID"])

    print(f"Combined dataset size: {len(combined_data)}")

    min_len = config["filtering"]["min_epitope_length"]
    max_len = config["filtering"]["max_epitope_length"]

    filtered_data = combined_data[
        combined_data["Modified_residues"].isna()
        & combined_data["Sequence"].str.len().between(min_len, max_len, inclusive="both")
    ].copy()

    print(f"Dataset size after basic filtering: {len(filtered_data)}")

    filtered_data = filtered_data.drop_duplicates(
        subset=["Protein_ID", "Sequence"]
    )

    print(f"Dataset size after removing duplicate seqs: {len(filtered_data)}")

    filtered_sequences = remove_nested_epitopes(filtered_data)

    print(f"Dataset size after removing nested epitopes: {len(filtered_sequences)}")

    output_path = DATA_DIR / "intermediate/wrangled_IEDB.csv"
    save_csv(filtered_sequences, output_path)

    print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()