#%%
import pandas as pd 
from Bio import SeqIO
import re
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
from pathlib import Path
import random

# Load combined FASTA
all_fasta = list(SeqIO.parse("../Data/all_proteomes.fasta", "fasta"))

def extract_genus_species_and_strain(organism_source):
    """
    Extracts genus + species as the first two words, and the rest as strain (if present).
    """
    cleaned = re.sub(r"\(.*?\)", "", organism_source).strip()
    parts = cleaned.split()

    if len(parts) >= 2:
        genus_species = ' '.join(parts[:2])
        strain = ' '.join(parts[2:]) if len(parts) > 2 else None
    else:
        genus_species = cleaned
        strain = "Unknown strain" 

    return genus_species, strain


def parse_fasta_to_df(fasta, dataset_name): 
    """
    Parse through FASTA records and extract metadata into a DataFrame.
    """
    metadata = []
    
    print(f"Processing {dataset_name}...")

    for seq_record in tqdm(fasta, desc=f"Parsing {dataset_name}", unit=" sequence"):
        header = seq_record.description

        # Extract the protein ID
        protein_id_match = re.search(r'\|([^|]+)\|', header)
        protein_id = protein_id_match.group(1) if protein_id_match else seq_record.id

        # Extract the full organism name
        organism_source = header.split('OS=')[1].split(' OX=')[0] if "OS=" in header else None

        # Extract genus/species and strain
        genus_species, strain = extract_genus_species_and_strain(organism_source) if organism_source else (None, "unknown strain")

        # Extract annotation
        annotation = None
        header_parts = header.split()
        if len(header_parts) > 1:
            possible_annotation = ' '.join(header_parts[1:])
            annotation = possible_annotation.split('OS=')[0].strip()

        # Extract gene name (GN=...)
        gene_name_match = re.search(r'GN=([^\s]+)', header)
        gene_name = gene_name_match.group(1) if gene_name_match else None


        # Shuffle protein sequence for null model
        seq_list = list(str(seq_record.seq))
        random.shuffle(seq_list)
        random_seq = str("".join(seq_list))

        # Add extracted data
        metadata.append([
            protein_id, genus_species, strain, annotation,
            gene_name, str(seq_record.seq), random_seq
        ])
    
    # Convert to DataFrame
    metadata_df = pd.DataFrame(metadata, columns=[
        "Protein_ID", "Genus_Species", "Strain", "Annotation",
        "pathogen_gene_name", "Sequence", "Random_sequence"
    ])
    
    return metadata_df


# Run the function
all_df = parse_fasta_to_df(all_fasta, "All Proteins")

# Save to CSV
output_csv_path = "../Data/wrangled_all_pathogen_prots.csv"
all_df.to_csv(output_csv_path, index=False)
print(f"Metadata saved to: {output_csv_path}")

# %%
