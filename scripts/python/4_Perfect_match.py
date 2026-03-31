#%% IMPORTS
import pandas as pd
import ahocorasick
from tqdm import tqdm
from helpers import DATA_DIR, load_csv
from collections import defaultdict

# load data
print("Loading data...")

pathogen_data_path = DATA_DIR / "intermediate/protein_metadata.csv"
IEDB_data_path = DATA_DIR / "intermediate/wrangled_IEDB.csv"
output_path = DATA_DIR / "proccesed/perfect_matches_2_0.csv"

pathogen_data = load_csv(pathogen_data_path)
pathogen_data = pathogen_data.rename(columns={
    "Genus_Species": "Organism",
})

IEDB_data = load_csv(IEDB_data_path)

# prepare epitopes
IEDB_epitopes = list(zip(
    IEDB_data["Assay_ID"],
    IEDB_data["Protein_source"],
    IEDB_data["Disease"],
    IEDB_data["Protein_ID"],
    IEDB_data["Sequence"],
    IEDB_data["epitope_start_pos"],
    IEDB_data["epitope_end_pos"]
))

# Build Aho-Corasick automaton for epitope matching
def build_automaton(epitopes):
    A = ahocorasick.Automaton()
    sub_map = defaultdict(list)

    for assay_id, source, disease, prot_id, epitope, start, end in epitopes:
        # generate all sub-epitopes from full length down to 9 aa
        for length in range(len(epitope), 8, -1):
            for i in range(len(epitope) - length + 1):
                sub = epitope[i:i+length]

                # Store metadata for each sub-epitope
                sub_map[sub].append((
                    assay_id, source, disease,
                    prot_id, sub, length,
                    start, end
                ))

    # Add sub-epitopes to automaton with their metadata into list, so duplicate sequences meta data doesnt get lost
    for sub, values in sub_map.items():
        A.add_word(sub, values)

    A.make_automaton()
    return A

print("Building automaton...")
A = build_automaton(IEDB_epitopes)

# Find matches between pathogen sequences and epitopes
def find_matches(pathogen_data, automaton):

    matches = []
    total_matches = 0

    print("Matching sequences...")

    for row in tqdm(pathogen_data.itertuples(index=False), total=len(pathogen_data)):

        seq = row.Sequence
        pid = row.Protein_ID
        proteome_id = getattr(row, "Proteome_ID", None)
        organism_id = getattr(row, "Organism_ID", None)
        pathogen_organism = row.Organism
        organism_source_fasta = getattr(row, "Organism_Source_FASTA", None)
        annot = row.Annotation
        gene = getattr(row, "Gene_name", None)

        for end_idx, hit_list in automaton.iter(seq):
            for (assay_id, source, disease,
                 ep_prot_id, sub, match_len,
                 ep_start, ep_end) in hit_list:

                match_start = end_idx - match_len + 2  # 1-based
                match_end = end_idx + 1

                matches.append((
                    assay_id, source, disease, ep_prot_id,
                    pid, proteome_id, organism_id,
                    pathogen_organism, organism_source_fasta,
                    annot, gene,
                    sub, match_len,
                    match_start, match_end,
                    ep_start, ep_end
                ))

                total_matches += 1

    match_df = pd.DataFrame(matches, columns=[
        "Assay_ID", "Epitope_Source", "Disease", "IEDB_Protein_ID",
        "Pathogen_Protein_ID", "Proteome_ID", "Organism_ID",
        "Pathogen_Organism", "Organism_Source_FASTA",
        "Pathogen_Annotation", "Pathogen_Gene_Name",
        "Matched_seq", "Match_Length",
        "Pathogen_Start", "Pathogen_End",
        "Epitope_Start", "Epitope_End"
    ])

    return match_df, total_matches


match_df, total_matches = find_matches(pathogen_data, A)

# keep only non-overlapping matches within each protein
def keep_non_overlapping_hits(group: pd.DataFrame) -> pd.DataFrame:
    """
    Keeps longest matches and removes overlapping shorter ones
    within each protein.
    """
    # Group values based first on match lenght, if equal then on start pos...
    group = group.sort_values(
        ["Match_Length", "Pathogen_Start", "Pathogen_End"],
        ascending=[False, True, True]
    )

    kept_rows = []
    kept_intervals = []

    # Itterate over matches in order of longest to shortest
    for _, row in group.iterrows():
        start = row["Pathogen_Start"]
        end = row["Pathogen_End"]

        # Check if current match overlaps with any already kept match
        overlaps = any(
            not (end < ks or start > ke)
            for ks, ke in kept_intervals
        )

        if not overlaps:
            kept_rows.append(row)
            kept_intervals.append((start, end))

    return pd.DataFrame(kept_rows)

print(f"Total matches before overlap filtering: {len(match_df)}")

grouped_results = []

# Group matches by these cols, to get relevant hits
group_cols = ["Assay_ID", "Pathogen_Protein_ID", "Proteome_ID", "Pathogen_Organism"]

for keys, group in match_df.groupby(group_cols, dropna=False, sort=False):
    filtered_group = keep_non_overlapping_hits(group.copy())

    # put grouping columns back 
    for col, val in zip(group_cols, keys if isinstance(keys, tuple) else (keys,)):
        filtered_group[col] = val

    grouped_results.append(filtered_group)

match_df = pd.concat(grouped_results, ignore_index=True)

print(f"Total matches after overlap filtering: {len(match_df)}")

# Random control: count matches in random sequences (should be much lower)
def count_random_matches(pathogen_data, automaton):
    count = 0

    for seq in tqdm(pathogen_data["Random_sequence"], desc="Random control"):
        for _ in automaton.iter(seq):
            count += 1
            break  # count presence only (fast)

    return count


random_matches = count_random_matches(pathogen_data, A)
print(f"Random sequence matches: {random_matches}")

# merge back with all epitopes to get full context and fill missing matches
all_epitopes = IEDB_data[[
    "Assay_ID", "Protein_source", "Disease",
    "Protein_ID", "Sequence",
    "epitope_start_pos", "epitope_end_pos"
]].rename(columns={
    "Protein_source": "Epitope_Source",
    "Protein_ID": "IEDB_Protein_ID"
}).drop_duplicates()

# Put back in the unmatched epitopes
full_result = all_epitopes.merge(
    match_df,
    how="left",
    on=["Assay_ID", "Epitope_Source", "Disease", "IEDB_Protein_ID"]
)

# Add a column to indicate if there was a match or not
full_result["Matched"] = ~full_result["Pathogen_Protein_ID"].isna()

# Only calculate end pos for matched rows
full_result.loc[full_result["Matched"], "Pathogen_End"] = (
    full_result.loc[full_result["Matched"], "Pathogen_Start"] +
    full_result.loc[full_result["Matched"], "Match_Length"] - 1
)


print(full_result["Matched"].value_counts(dropna=False))
# --------------------------------------------

# save results
full_result.to_csv(output_path, index=False)
print(f"Saved results to: {output_path}")
# %%