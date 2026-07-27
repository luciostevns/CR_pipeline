import re

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

from crossreactivity.io import DATA_DIR, load_csv, save_csv


def load_input_data():

    fasta_records = list(
        SeqIO.parse(DATA_DIR / "raw/all_proteomes.fasta", "fasta")
    )

    pathogenic_df = load_csv(
        DATA_DIR / "intermediate/pathogenic_bacteria_proteome.csv"
    )

    return fasta_records, pathogenic_df

def parse_header_annotation(header: str) -> str | None:
    """
    Extract a real protein annotation from a FASTA header.

    For UniProtKB:
        sp|P0A6F5|CH60_ECOLI Chaperonin GroEL OS=...
        -> Chaperonin GroEL

    For UniParc:
        UPI00000BD9A2 status=active PROTEOME=... SOURCE_DB=UniParc
        -> None
    """

    # Remove the first token, which is the sequence ID
    parts = header.split(maxsplit=1)

    if len(parts) == 1:
        return None

    rest = parts[1]

    # Cut away known metadata fields
    for marker in [" OS=", " OX=", " GN=", " PE=", " SV=", " PROTEOME=", " SOURCE_DB=", " status="]:
        if marker in rest:
            rest = rest.split(marker)[0].strip()

    if not rest:
        return None

    # Avoid storing metadata as an annotation
    bad_prefixes = (
        "status=",
        "PROTEOME=",
        "SOURCE_DB=",
    )

    if rest.startswith(bad_prefixes):
        return None

    return rest

def parse_fasta_to_df(fasta_records, dataset_name: str) -> pd.DataFrame:
    """
    Parse FASTA records and extract protein-level metadata from FASTA headers.
    """
    metadata = []

    for seq_record in tqdm(
        fasta_records,
        desc=f"Parsing {dataset_name}",
        unit="seq"
    ):
        header = seq_record.description

        # Protein ID
        protein_id_match = re.search(r"\|([^|]+)\|", header)
        protein_id = protein_id_match.group(1) if protein_id_match else seq_record.id

        # Proteome ID
        proteome_match = re.search(r"PROTEOME=([^\s]+)", header)
        proteome_id = proteome_match.group(1) if proteome_match else None

        # Protein annotation
        annotation = parse_header_annotation(header)

        # Gene name
        gene_name_match = re.search(r"GN=([^\s]+)", header)
        gene_name = gene_name_match.group(1) if gene_name_match else None

        metadata.append(
            [
                protein_id,
                proteome_id,
                annotation,
                gene_name,
                str(seq_record.seq),
            ]
        )

    metadata_df = pd.DataFrame(
        metadata,
        columns=[
            "Protein_ID",
            "Proteome_ID",
            "Annotation",
            "Gene_name",
            "Sequence",
        ]
    )

    return metadata_df


def merge_proteome_metadata(protein_df: pd.DataFrame, pathogenic_df: pd.DataFrame) -> pd.DataFrame:

    protein_metadata_df = protein_df.merge(
        pathogenic_df[
            [
                "Proteome_ID",
                "Scientific_name",
                "Genus_species",
                "Strain",
            ]
        ],
        how="left",
        on="Proteome_ID"
    )

    return protein_metadata_df


def filter_rows_without_strain(protein_metadata_df: pd.DataFrame) -> pd.DataFrame:
    before = len(protein_metadata_df)

    protein_metadata_df = protein_metadata_df[
        protein_metadata_df["Strain"].notna()
        & (protein_metadata_df["Strain"].str.strip() != "")
    ].copy()

    after = len(protein_metadata_df)

    print(
        f"Removed {before - after} rows without strain info "
        f"({(before - after) / before:.2%})"
    )

    return protein_metadata_df

def remove_duplicate_proteins(
    protein_metadata_df: pd.DataFrame,
) -> pd.DataFrame:
    """Remove duplicate proteins and identify occurrences across species/genera."""
    print("Checking for duplicate proteins...")

    protein_metadata_df = protein_metadata_df.copy()

    keys = pd.DataFrame(
        {
            "protein": pd.factorize(
                protein_metadata_df["Protein_ID"]
            )[0],
            "sequence": pd.factorize(
                protein_metadata_df["Sequence"]
            )[0],
            "species": pd.factorize(
                protein_metadata_df["Genus_species"]
            )[0],
        },
        index=protein_metadata_df.index,
    )

    protein_key = ["protein", "sequence"]
    duplicate_key = protein_key + ["species"]

    # Identify and remove repeated rows within the same species
    duplicated = keys.duplicated(duplicate_key, keep=False)
    keep_mask = ~keys.duplicated(duplicate_key, keep="first")

    result = protein_metadata_df.loc[keep_mask].copy()
    retained_keys = keys.loc[keep_mask].copy()

    retained_keys["genus"] = pd.factorize(
        result["Genus_species"].str.split().str[0]
    )[0]

    groups = retained_keys.groupby(protein_key, sort=False)

    # Same Protein_ID + Sequence occurring in multiple species/genera
    cross_species_mask = groups["species"].transform("size").gt(1)
    cross_genus_mask = groups["genus"].transform("nunique").gt(1)
    first_occurrence = ~retained_keys.duplicated(protein_key)

    # Repeated rows found only within one species
    same_species_mask = (
        duplicated.loc[keep_mask].to_numpy()
        & ~cross_species_mask.to_numpy()
    )

    result["Duplicate_across_genus_species"] = (
        cross_species_mask.to_numpy()
    )
    result["Duplicate_across_genera"] = (
        cross_genus_mask.to_numpy()
    )

    print(
        f"Removed duplicate protein rows: "
        f"{len(protein_metadata_df) - len(result):,}"
    )
    print(
        f"Unique combinations duplicated within one Genus_species: "
        f"{same_species_mask.sum():,}"
    )
    print(
        f"Unique combinations occurring across Genus_species: "
        f"{(cross_species_mask & first_occurrence).sum():,}"
    )
    print(
        f"Unique combinations occurring across genera: "
        f"{(cross_genus_mask & first_occurrence).sum():,}"
    )

    columns = ["Protein_ID", "Sequence", "Genus_species"]

    print("\nFive within-species examples:")
    print(
        result.loc[same_species_mask, columns]
        .assign(Sequence=lambda x: x["Sequence"].str[:30] + "...")
        .head(5)
        .to_string(index=False)
    )

    example_keys = retained_keys.loc[
        cross_species_mask & first_occurrence,
        protein_key,
    ].head(5)

    selected = pd.MultiIndex.from_frame(example_keys)
    all_retained_keys = pd.MultiIndex.from_frame(
        retained_keys[protein_key]
    )

    cross_examples = (
        result.loc[all_retained_keys.isin(selected), columns]
        .groupby(
            ["Protein_ID", "Sequence"],
            sort=False,
            dropna=False,
        )
        .agg(
            Genus_species=(
                "Genus_species",
                lambda x: ", ".join(x.dropna().unique()),
            )
        )
        .reset_index()
        .assign(Sequence=lambda x: x["Sequence"].str[:30] + "...")
    )

    print("\nFive cross-species examples:")
    print(cross_examples.to_string(index=False))

    genus_example_keys = retained_keys.loc[
        cross_genus_mask & first_occurrence,
        protein_key,
    ].head(5)

    selected = pd.MultiIndex.from_frame(genus_example_keys)

    cross_genus_examples = (
        result.loc[all_retained_keys.isin(selected), columns]
        .groupby(
            ["Protein_ID", "Sequence"],
            sort=False,
            dropna=False,
        )
        .agg(
            n_species=("Genus_species", "nunique"),
            genera=(
                "Genus_species",
                lambda x: ", ".join(
                    sorted(
                        x.dropna()
                        .str.split()
                        .str[0]
                        .unique()
                    )
                ),
            ),
        )
        .reset_index()
        .assign(Sequence=lambda x: x["Sequence"].str[:30] + "...")
    )

    print("\nFive cross-genus examples:")
    print(cross_genus_examples.to_string(index=False))

    # Get the Protein_ID + Sequence combinations occurring
    # across species and across genera.
    cross_species_keys = (
        retained_keys.loc[cross_species_mask, protein_key]
        .drop_duplicates()
    )

    cross_genus_keys = (
        retained_keys.loc[cross_genus_mask, protein_key]
        .drop_duplicates()
    )

    # Apply these keys to the original dataframe so annotations
    # are also collected from rows removed as duplicates.
    original_keys = pd.MultiIndex.from_frame(
        keys[protein_key]
    )
    species_keys = pd.MultiIndex.from_frame(
        cross_species_keys
    )
    genus_keys = pd.MultiIndex.from_frame(
        cross_genus_keys
    )

    original_cross_species_mask = original_keys.isin(
        species_keys
    )
    original_cross_genus_mask = original_keys.isin(
        genus_keys
    )

    species_rows = protein_metadata_df.loc[original_cross_species_mask].copy()
    genus_rows = protein_metadata_df.loc[original_cross_genus_mask].copy()


    def check_annotations(rows: pd.DataFrame, label: str) -> None:
        cleaned_annotations = (
            rows["Annotation"]
            .fillna("")
            .astype(str)
            .str.strip()
        )

        print(f"\nAnnotation check for {label}:")
        print(f"Total rows: {len(rows):,}")
        print(f"Annotation is NA: {rows['Annotation'].isna().sum():,}")
        print(f"Annotation is blank: {cleaned_annotations.eq('').sum():,}")
        print(f"Annotation is available: {cleaned_annotations.ne('').sum():,}")

        print(
            rows[
                [
                    "Protein_ID",
                    "Proteome_ID",
                    "Genus_species",
                    "Annotation",
                ]
            ]
            .head(10)
            .to_string(index=False)
        )


    check_annotations(species_rows, "cross-species proteins")
    check_annotations(genus_rows, "cross-genus proteins")

    duplicate_annotations_species = (
        protein_metadata_df.loc[
            original_cross_species_mask,
            ["Annotation"],
        ]
        .dropna(subset=["Annotation"])
        .assign(
            Annotation=lambda df: df["Annotation"].str.strip()
        )
        .query("Annotation != ''")
        .drop_duplicates()
        .sort_values("Annotation")
        .reset_index(drop=True)
    )

    duplicate_annotations_genera = (
        protein_metadata_df.loc[
            original_cross_genus_mask,
            ["Annotation"],
        ]
        .dropna(subset=["Annotation"])
        .assign(
            Annotation=lambda df: df["Annotation"].str.strip()
        )
        .query("Annotation != ''")
        .drop_duplicates()
        .sort_values("Annotation")
        .reset_index(drop=True)
    )

    save_csv(
        duplicate_annotations_species,
        DATA_DIR
        / "intermediate/duplicate_annotations_across_species.csv",
    )

    save_csv(
        duplicate_annotations_genera,
        DATA_DIR
        / "intermediate/duplicate_annotations_across_genera.csv",
    )

    print(
        f"Original rows belonging to cross-species proteins: "
        f"{original_cross_species_mask.sum():,}"
    )
    print(
        f"Unique annotations occurring across species: "
        f"{len(duplicate_annotations_species):,}"
    )
    print(
        f"Original rows belonging to cross-genus proteins: "
        f"{original_cross_genus_mask.sum():,}"
    )
    print(
        f"Unique annotations occurring across genera: "
        f"{len(duplicate_annotations_genera):,}"
    )

    return result


def main() -> None:
    fasta_records, pathogenic_df = load_input_data()

    protein_df = parse_fasta_to_df(
        fasta_records=fasta_records,
        dataset_name="All Proteins"
    )

    protein_metadata_df = merge_proteome_metadata(
        protein_df=protein_df,
        pathogenic_df=pathogenic_df
    )

    protein_metadata_df = filter_rows_without_strain(
        protein_metadata_df=protein_metadata_df
    )

    protein_metadata_df = remove_duplicate_proteins(
        protein_metadata_df=protein_metadata_df
    )

    output_path = DATA_DIR / "intermediate/protein_metadata.csv"

    save_csv(protein_metadata_df, output_path)


if __name__ == "__main__":
    main()