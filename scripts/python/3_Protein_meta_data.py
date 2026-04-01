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

        # Randomized sequence
        seq_list = list(str(seq_record.seq))
        random.shuffle(seq_list)
        random_seq = "".join(seq_list)

        metadata.append([
            protein_id,
            proteome_id,
            annotation,
            gene_name,
            str(seq_record.seq),
            random_seq
        ])

    metadata_df = pd.DataFrame(metadata, columns=[
        "Protein_ID",
        "Proteome_ID",
        "Annotation",
        "Gene_name",
        "Sequence",
        "Random_sequence"
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

# Diagnostics
print("\nMissing values after merge:")
print(merged_df[["Scientific_name", "Genus_species", "Strain"]].isna().sum())

print("\nUnique counts:")
print(f"Unique Proteome_ID in FASTA: {merged_df['Proteome_ID'].nunique()}")
print(f"Unique Genus_species: {merged_df['Genus_species'].nunique()}")
print(f"Unique Strain: {merged_df['Strain'].nunique()}")

# Save result
print("Saving metadata to CSV...")
save_csv(merged_df, DATA_DIR / "intermediate/protein_metadata.csv")

print("Done")