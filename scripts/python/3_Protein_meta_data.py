import pandas as pd
from Bio import SeqIO
import re
from tqdm import tqdm
import random
from helpers import DATA_DIR, save_csv, load_csv

# Load data
print("Loading combined FASTA file...")
all_fasta = list(SeqIO.parse(DATA_DIR / "raw/all_proteomes.fasta", "fasta"))

print("Loading pathogenic proteome metadata...")
pathogenic_df = load_csv(DATA_DIR / "intermediate/pathogenic_bacteria_proteome.csv")


def parse_fasta_to_df(fasta, dataset_name):
    """
    Parse FASTA records and extract protein-level metadata from FASTA headers.
    """
    metadata = []

    print(f"Processing {dataset_name}...")

    for seq_record in tqdm(fasta, desc=f"Parsing {dataset_name}", unit="seq"):
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

        metadata.append([
            protein_id,
            proteome_id,
            annotation,
            gene_name,
            str(seq_record.seq)
        ])

    metadata_df = pd.DataFrame(metadata, columns=[
        "Protein_ID",
        "Proteome_ID",
        "Annotation",
        "Gene_name",
        "Sequence"
    ])

    return metadata_df


# Parse FASTA
all_df = parse_fasta_to_df(all_fasta, "All Proteins")

# Merge in proteome-level metadata
print("Merging in proteome-level metadata...")
merged_df = all_df.merge(
    pathogenic_df[
        [
            "Proteome_ID",
            "Scientific_name",
            "Genus_species",
            "Strain"
        ]
    ],
    how="left",
    on="Proteome_ID"
)

# Remove rows without strain information
before = len(merged_df)
merged_df = merged_df[
    merged_df["Strain"].notna() & (merged_df["Strain"].str.strip() != "")
].copy()
after = len(merged_df)

print(f"Removed {before - after} rows without strain info ({(before - after)/before:.2%})")

# Diagnostics
print("\nMissing values after merge:")
print(merged_df[["Scientific_name", "Genus_species", "Strain"]].isna().sum())

print(f"\nTotal proteins(rows) {after}")
print(f"\nTotal unique rows:\n{merged_df.nunique()}")

n_unique_rows = merged_df.drop_duplicates().shape[0]
print(f"\nUnique rows (all columns): {n_unique_rows}")

dup_counts = (
    merged_df.groupby("Sequence")
    .agg(
        count=("Sequence", "size"),
        annotations=("Annotation", lambda x: list(pd.unique(x)))
    )
    .sort_values("count", ascending=False)
    .reset_index()
)

print(dup_counts[["annotations", "count"]].head(10))

# Save result
print("Saving metadata to CSV...")
save_csv(merged_df, DATA_DIR / "intermediate/protein_metadata.csv")

print("Done")