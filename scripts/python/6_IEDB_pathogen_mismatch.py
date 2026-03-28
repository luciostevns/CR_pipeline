#%%
## IMPORTANT, RUN SCRIPT ONLY IF NEEDED SINCE EDLIB IS NOT DOWNLOADED

import pandas as pd
import numpy as np
from tqdm import tqdm
from collections import defaultdict
import edlib

# Load data
pathogen_data = pd.read_csv("../Data/wrangled_rep_pathogen_prots.csv")
IEDB_data = pd.read_csv("../Data/wrangled_IEDB.csv")

# Extract relevant columns
epitope_sequences = IEDB_data["Sequence"].to_numpy()
epitope_ids = IEDB_data["Assay_ID"].to_numpy()
epitope_sources = IEDB_data["Epitope - Molecule Parent"].to_numpy()

protein_sequences = pathogen_data["Sequence"].to_numpy()
protein_ids = pathogen_data["Protein ID"].to_numpy()
organism = pathogen_data["Organism Source"].to_numpy()

# Function to generate all 9-mers from a sequence (optimized)
def generate_9mers(seq):
    return [seq[i:i+9] for i in range(len(seq) - 8)] if len(seq) >= 9 else []

# Using edlib to calculate Hamming distance (by leveraging edit distance with a gap cost of -1 for insertions and deletions)
def hamming_distance(s1, s2, max_mismatches):
    alignment = edlib.align(s1, s2, mode="NW")
    mismatches = alignment['editDistance']
    return mismatches if mismatches <= max_mismatches else float('inf')

# Define allowed mismatches
max_mismatches = 4
prefix_length = max_mismatches  # Ensure correctness of comparisons

# Step 2: Group epitope 9-mers by prefix using defaultdict
epitope_dict = defaultdict(list)

for assay_id, epitope_source, seq in zip(epitope_ids, epitope_sources, epitope_sequences):
    for nine_mer in generate_9mers(seq):
        prefix = nine_mer[:prefix_length]
        epitope_dict[prefix].append((nine_mer, assay_id, epitope_source))

# Function to find matches in a protein sequence
def find_fuzzy_matches(protein, protein_id, organism, epitope_dict, max_mismatches):
    local_matches = []
    protein_9mers = generate_9mers(protein)
    
    for protein_9mer in protein_9mers:
        prefix = protein_9mer[:prefix_length]
        if prefix in epitope_dict:
            for epitope_9mer, assay_id, epitope_source in epitope_dict[prefix]:
                # Perform fuzzy matching only if prefix matches
                mismatches = hamming_distance(protein_9mer, epitope_9mer, max_mismatches)
                if mismatches <= max_mismatches:
                    local_matches.append((assay_id, epitope_source, protein_id, organism, protein_9mer, epitope_9mer))
    return local_matches

# Store results
matches = []

# Iterate over protein sequences and store matches
with tqdm(total=len(protein_sequences), desc="Processing Proteins") as pbar:
    for protein, protein_id, org in zip(protein_sequences, protein_ids, organism):
        matches.extend(find_fuzzy_matches(protein, protein_id, org, epitope_dict, max_mismatches))
        pbar.update(1)

# Convert to DataFrame for easier viewing & saving
match_df = pd.DataFrame(matches, columns=["Assay_ID", "Epitope Source", "Protein_ID", "Organism Source", "Matched_9mer", "Epitope_9mer"])

# Print total matches and head of DataFrame
print(f"\nTotal Matches Found: {len(matches)}")

# Optionally save to CSV
match_df.to_csv("../Data/matching_9mers.csv", index=False)

# %%
