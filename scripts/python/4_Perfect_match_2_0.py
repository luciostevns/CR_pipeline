#%% IMPORTS
import pandas as pd
import numpy as np
import ahocorasick
from tqdm import tqdm

#%% LOAD DATA
pathogen_data_path = "../Data/wrangled_all_pathogen_prots.csv"
IEDB_data_path = "../Data/wrangled_IEDB.csv"
output_path = "../Data/perfect_matches_2_0.csv"

# Load data
pathogen_data = pd.read_csv(pathogen_data_path)
IEDB_data = pd.read_csv(IEDB_data_path)

IEDB_epitopes = list(zip(
    IEDB_data["Assay_ID"],
    IEDB_data["Protein_source"],
    IEDB_data["Disease"],
    IEDB_data["Protein_ID"],
    IEDB_data["Sequence"],
    IEDB_data["epitope_start_pos"],
    IEDB_data["epitope_end_pos"]
))

#%% BUILD AHO-CORASICK AUTOMATON
def build_aho_corasick_automaton(epitopes):
    """
    Builds an Aho-Corasick automaton for efficient multi-pattern matching with a sliding window.
    """

    A = ahocorasick.Automaton()
    for assay_id, epitope_source, disease, epitope_protein_id, epitope, epitope_start, epitope_end in epitopes:
        for length in range(len(epitope), 8, -1):  # 15 → 9-mers
            for i in range(len(epitope) - length + 1):
                sub_epitope = epitope[i:i + length]
                A.add_word(sub_epitope, (
                    assay_id, epitope_source, disease,
                    epitope_protein_id, sub_epitope, length,
                    epitope_start, epitope_end
                ))
    A.make_automaton()
    return A

A = build_aho_corasick_automaton(IEDB_epitopes)

#%% FIND MATCHES
def find_matches(pathogen_data, automaton):
    """
    Finds matches between pathogen protein sequences and epitopes using the Aho-Corasick automaton.
    Returns a DataFrame with match details and metadata, and the total raw count of matches.
    """

    matches = []
    total_matches = 0

    seqs = pathogen_data["Sequence"].to_numpy()
    random_seqs = pathogen_data["Random_sequence"].to_numpy()
    pids = pathogen_data["Protein_ID"].to_numpy()
    orgs = pathogen_data["Genus_Species"].to_numpy()
    annots = pathogen_data["Annotation"].to_numpy()
    strains = pathogen_data["Strain"].to_numpy()
    genes = pathogen_data["pathogen_gene_name"].to_numpy()

    with tqdm(total=len(seqs), desc="Matching epitopes") as pbar:
        for seq, pid, org, strain, annot, gene in zip(seqs, pids, orgs, strains, annots, genes):
            row_matches = []

            for end_idx, (assay_id, epitope_source, disease,
                          epitope_protein_id, sub_epitope, match_len,
                          epitope_start, epitope_end) in automaton.iter(seq):
                match_start = end_idx - match_len + 2  # 1-based
                match_end = end_idx + 1

                row_matches.append((assay_id, epitope_source, disease, epitope_protein_id,
                                    pid, org, strain, annot, gene,
                                    sub_epitope, match_len,
                                    match_start, match_end,
                                    epitope_start, epitope_end))
                total_matches += 1

            if row_matches:
                row_df = pd.DataFrame(row_matches, columns=[
                    "Assay_ID", "Epitope_Source", "Disease", "IEDB_Protein_ID",
                    "Pathogen_Protein_ID", "Organism_Source", "Strain", "Pathogen_Annotation",
                    "Pathogen_Gene_Name", "Matched_9mer", "Match_Length",
                    "Pathogen_Protein_Start_Pos", "Pathogen_Protein_End_Pos",
                    "Epitope_Start_Pos", "Epitope_End_Pos"
                ])
                row_df = row_df.sort_values("Match_Length", ascending=False)
                row_df = row_df.drop_duplicates(subset=["Assay_ID", "Pathogen_Protein_ID", "Strain"])
                matches.extend(row_df.to_records(index=False))

            pbar.update(1)

    match_df = pd.DataFrame.from_records(matches, columns=[
        "Assay_ID", "Epitope_Source", "Disease", "IEDB_Protein_ID",
        "Pathogen_Protein_ID", "Organism_Source", "Strain", "Pathogen_Annotation",
        "Pathogen_Gene_Name", "Matched_9mer", "Match_Length",
        "Pathogen_Protein_Start_Pos", "Pathogen_Protein_End_Pos",
        "Epitope_Start_Pos", "Epitope_End_Pos"
    ])

    match_df = match_df.sort_values("Match_Length", ascending=False)
    match_df = match_df.drop_duplicates(subset=["Assay_ID", "Pathogen_Protein_ID", "Strain"])

    return match_df, total_matches

#%% EXECUTE MATCHING
match_df, total_matches = find_matches(pathogen_data, A)
print(f"Total raw matches in sequences: {total_matches}")

#%% COUNT MATCHES IN RANDOM SEQUENCES
def count_matches_in_random_seqs(pathogen_data, automaton):
    """
    Counts the total number of matches found in random sequences using the Aho-Corasick automaton.
    """
    count = 0
    random_seqs = pathogen_data["Random_sequence"].to_numpy()
    for seq in random_seqs:
        for end_idx, match_data in automaton.iter(seq):
            count += 1
    return count

random_matches_count = count_matches_in_random_seqs(pathogen_data, A)
print(f"Total matches in random sequences: {random_matches_count}")

#%% MERGE WITH ALL EPITOPES (even unmatched ones)
all_epitopes = IEDB_data[[
    "Assay_ID", "Protein_source", "Disease", "Protein_ID", "Sequence",
    "epitope_start_pos", "epitope_end_pos"
]].rename(columns={
    "Protein_source": "Epitope_Source",
    "Protein_ID": "IEDB_Protein_ID"
}).drop_duplicates()

full_result = all_epitopes.merge(
    match_df,
    how="left",
    on=["Assay_ID", "Epitope_Source", "Disease", "IEDB_Protein_ID"]
)


full_result["Match_Length"] = full_result["Match_Length"].fillna(full_result["Sequence"].str.len())
full_result["Pathogen_Protein_End_Pos"] = full_result["Pathogen_Protein_Start_Pos"] + full_result["Match_Length"] - 1

full_result["Matched"] = ~full_result["Pathogen_Protein_ID"].isna()

# Save final result
full_result.to_csv(output_path, index=False)
print(f"✅ Saved full epitope match table (including unmatched) to: {output_path}")
# %%

