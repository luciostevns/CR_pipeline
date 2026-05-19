import time

import pandas as pd
import requests
from tqdm import tqdm

from crossreactivity.io import (
    DATA_DIR,
    load_tsv,
    load_tsv_robust,
    save_csv
)
from crossreactivity.uniprot import fetch_proteome_metadata


def load_input_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    pathogen_ref = load_tsv_robust(DATA_DIR / "raw/isolates.tsv")

    referenced_proteome = load_tsv(
        DATA_DIR / "raw/Uniprot_raw_refferenced_bacterial_proteomes.tsv"
    )

    other_proteome = load_tsv(
        DATA_DIR / "raw/Uniprot_raw_unrefferenced_bacterial_proteomes.tsv"
    )

    proteome = pd.concat(
        [referenced_proteome, other_proteome],
        ignore_index=True
    )

    return pathogen_ref, proteome


def clean_merge_ids(
    pathogen_ref: pd.DataFrame,
    proteome: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:

    proteome = proteome.dropna(
        subset=["Proteome Id", "Genome assembly ID"]
    ).copy()

    pathogen_ref = pathogen_ref.dropna(
        subset=["Assembly"]
    ).copy()

    proteome["Genome assembly ID"] = (
        proteome["Genome assembly ID"]
        .astype(str)
        .str.strip()
        .str.upper()
    )

    pathogen_ref["Assembly"] = (
        pathogen_ref["Assembly"]
        .astype(str)
        .str.strip()
        .str.upper()
    )

    return pathogen_ref, proteome


def filter_pathogenic_proteomes(
    pathogen_ref: pd.DataFrame,
    proteome: pd.DataFrame
) -> pd.DataFrame:

    valid_assemblies = set(pathogen_ref["Assembly"])

    pathogenic_proteomes = proteome[
        proteome["Genome assembly ID"].isin(valid_assemblies)
    ].copy()

    return pathogenic_proteomes


def make_base_proteome_metadata(
    pathogenic_proteomes: pd.DataFrame
) -> pd.DataFrame:

    proteome_metadata = pathogenic_proteomes[
        ["Proteome Id", "Genome assembly ID", "Protein count"]
    ].copy()

    proteome_metadata = proteome_metadata.rename(
        columns={
            "Proteome Id": "Proteome_ID",
            "Genome assembly ID": "Genome_Assembly_ID",
            "Protein count": "Protein_Count",
        }
    )

    return proteome_metadata


def fetch_all_proteome_metadata(
    proteome_ids: pd.Series
) -> pd.DataFrame:

    rest_metadata = []
    session = requests.Session()

    print("Fetching Scientific_name / Genus_species / Strain from UniProt REST...")

    for proteome_id in tqdm(
        proteome_ids,
        desc="Fetching proteome metadata"
    ):

        meta = fetch_proteome_metadata(
            proteome_id,
            session=session
        )

        rest_metadata.append(meta)

        time.sleep(0.2)

    return pd.DataFrame(rest_metadata)


def main() -> None:

    pathogen_ref, proteome = load_input_data()

    pathogen_ref, proteome = clean_merge_ids(
        pathogen_ref=pathogen_ref,
        proteome=proteome
    )

    print(f"Total proteome entries before filtering: {len(proteome)}")

    pathogenic_proteomes = filter_pathogenic_proteomes(
        pathogen_ref=pathogen_ref,
        proteome=proteome
    )

    print(
        "Proteome entries after filtering for pathogenic assemblies: "
        f"{len(pathogenic_proteomes)}"
    )

    proteome_metadata = make_base_proteome_metadata(
        pathogenic_proteomes
    )

    rest_metadata_df = fetch_all_proteome_metadata(
        proteome_metadata["Proteome_ID"]
    )

    proteome_metadata = proteome_metadata.merge(
        rest_metadata_df,
        how="left",
        on="Proteome_ID"
    )

    proteome_metadata = proteome_metadata[
        [
            "Proteome_ID",
            "Genome_Assembly_ID",
            "Protein_Count",
            "Scientific_name",
            "Genus_species",
            "Strain",
        ]
    ]

    pathogenic_output = (
        DATA_DIR / "intermediate/pathogenic_bacteria_proteome.csv"
    )

    proteome_ids_output = (
        DATA_DIR / "intermediate/proteome_ids.csv"
    )

    save_csv(
        proteome_metadata,
        pathogenic_output
    )

    save_csv(
        proteome_metadata[["Proteome_ID"]],
        proteome_ids_output
    )

    total_protein_count = pd.to_numeric(
        proteome_metadata["Protein_Count"],
        errors="coerce"
    ).sum()

    matching_ids = (
        set(proteome["Genome assembly ID"])
        & set(pathogen_ref["Assembly"])
    )

    print(f"Total annotated protein count: {total_protein_count}")
    print(f"Number of matching Assembly IDs: {len(matching_ids)}")

    print("\nMissing values:")
    print(
        proteome_metadata[
            ["Scientific_name", "Genus_species", "Strain"]
        ].isna().sum()
    )

    print("Saved:")
    print(f"- {pathogenic_output}")
    print(f"- {proteome_ids_output}")


if __name__ == "__main__":
    main()