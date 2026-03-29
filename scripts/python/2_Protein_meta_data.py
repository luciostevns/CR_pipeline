import pandas as pd
from Bio import SeqIO
import re
from tqdm import tqdm
import random
from helpers import DATA_DIR, save_csv, load_csv

# Load data
print("Loading combined FASTA file...")
all_fasta = list(SeqIO.parse(DATA_DIR / "raw/all_proteomes.fasta", "fasta"))

print("Loading pathogenic bacteria metadata...")
pathogenic_df = load_csv(DATA_DIR / "intermediate/pathogenic_bacteria_proteome.csv")


# FASTA parsing function
def parse_fasta_to_df(fasta, dataset_name):
    """
    Parse FASTA records and extract metadata.
    """
    metadata = []

    print(f"Processing {dataset_name}...")

    for seq_record in tqdm(fasta, desc=f"Parsing {dataset_name}", unit="seq"):
        header = seq_record.description

        # Protein ID
        protein_id_match = re.search(r'\|([^|]+)\|', header)
        protein_id = protein_id_match.group(1) if protein_id_match else seq_record.id

        # Proteome ID
        proteome_match = re.search(r'PROTEOME=([^\s]+)', header)
        proteome_id = proteome_match.group(1) if proteome_match else None

        # protein annotation
        annotation = None
        header_parts = header.split()
        if len(header_parts) > 1:
            annotation = ' '.join(header_parts[1:]).split('OS=')[0].strip()

        # Gene name
        gene_name_match = re.search(r'GN=([^\s]+)', header)
        gene_name = gene_name_match.group(1) if gene_name_match else None

        # Randomized sequence (null model)
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


# Run parsing
all_df = parse_fasta_to_df(all_fasta, "All Proteins")


# Merge with metadata
print("Merging with pathogenic metadata...")

merged_df = all_df.merge(
    pathogenic_df,
    how="left",
    left_on="Proteome_ID",
    right_on="Proteome Id"
)

cols_to_drop = [
    "Assembly",
    "Protein count"
]

merged_df = merged_df.drop(columns=cols_to_drop, errors="ignore")

# Diagnostics
missing = merged_df["Proteome Id"].isna().sum()
print(f"Missing Proteome_ID matches: {missing}")


# Save result
print("Saving metadata to CSV...")
save_csv(merged_df, DATA_DIR / "intermediate/protein_metadata.csv")

print("Done")