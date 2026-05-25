import re

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

from crossreactivity.io import DATA_DIR, load_csv, save_csv


def load_input_data():
    print("Loading combined FASTA file...")
    fasta_records = list(
        SeqIO.parse(DATA_DIR / "raw/all_proteomes.fasta", "fasta")
    )

    print("Loading pathogenic proteome metadata...")
    pathogenic_df = load_csv(
        DATA_DIR / "intermediate/pathogenic_bacteria_proteome.csv"
    )

    return fasta_records, pathogenic_df


def parse_fasta_to_df(fasta_records, dataset_name: str) -> pd.DataFrame:
    """
    Parse FASTA records and extract protein-level metadata from FASTA headers.
    """
    metadata = []

    print(f"Processing {dataset_name}...")

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
        annotation = None
        header_parts = header.split()

        if len(header_parts) > 1:
            annotation = " ".join(header_parts[1:]).split("OS=")[0].strip()

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
    print("Merging in proteome-level metadata...")

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


def print_diagnostics(protein_metadata_df: pd.DataFrame) -> None:
    print("\nMissing values after merge:")
    print(
        protein_metadata_df[
            ["Scientific_name", "Genus_species", "Strain"]
        ].isna().sum()
    )

    print(f"\nTotal proteins(rows): {len(protein_metadata_df)}")
    print(f"\nTotal unique values per column:\n{protein_metadata_df.nunique()}")

    n_unique_rows = protein_metadata_df.drop_duplicates().shape[0]
    print(f"\nUnique rows (all columns): {n_unique_rows}")

    duplicate_sequence_counts = (
        protein_metadata_df
        .groupby("Sequence")
        .agg(
            count=("Sequence", "size"),
            annotations=("Annotation", lambda x: list(pd.unique(x)))
        )
        .sort_values("count", ascending=False)
        .reset_index()
    )

    print("\nMost duplicated protein sequences:")
    print(duplicate_sequence_counts[["annotations", "count"]].head(10))


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

    print_diagnostics(protein_metadata_df)

    output_path = DATA_DIR / "intermediate/protein_metadata.csv"

    print("Saving metadata to CSV...")
    save_csv(protein_metadata_df, output_path)

    print(f"Saved: {output_path}")
    print("Done")


if __name__ == "__main__":
    main()