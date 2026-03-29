#%% IMPORTS
import pandas as pd
import ahocorasick
from tqdm import tqdm
from helpers import DATA_DIR, load_csv

# load data
print("Loading data...")

pathogen_data_path = DATA_DIR / "intermediate/protein_metadata.csv"
IEDB_data_path = DATA_DIR / "intermediate/wrangled_IEDB.csv"
output_path = DATA_DIR / "proccesed/perfect_matches_2_0.csv"

pathogen_data = load_csv(pathogen_data_path)
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

    for assay_id, source, disease, prot_id, epitope, start, end in epitopes:
        for length in range(len(epitope), 8, -1):  # 15 → 9-mers
            for i in range(len(epitope) - length + 1):
                sub = epitope[i:i+length]

                A.add_word(sub, (
                    assay_id, source, disease,
                    prot_id, sub, length,
                    start, end
                ))

    A.make_automaton()
    return A

print("Building automaton...")
A = build_automaton(IEDB_epitopes)
print(f"Automaton size: {len(A)}")

# Find matches between pathogen sequences and epitopes
def find_matches(pathogen_data, automaton):

    matches = []
    total_matches = 0

    print("Matching sequences...")

    for row in tqdm(pathogen_data.itertuples(index=False), total=len(pathogen_data)):

        seq = row.Sequence
        pid = row.Protein_ID
        strain = row.Strain
        annot = row.Annotation
        gene = getattr(row, "Gene_name", None)

        for end_idx, (assay_id, source, disease,
                      ep_prot_id, sub, match_len,
                      ep_start, ep_end) in automaton.iter(seq):

            match_start = end_idx - match_len + 2  # 1-based
            match_end = end_idx + 1

            matches.append((
                assay_id, source, disease, ep_prot_id,
                pid, strain, annot, gene,
                sub, match_len,
                match_start, match_end,
                ep_start, ep_end
            ))

            total_matches += 1

    match_df = pd.DataFrame(matches, columns=[
        "Assay_ID", "Epitope_Source", "Disease", "IEDB_Protein_ID",
        "Pathogen_Protein_ID", "Strain", "Pathogen_Annotation",
        "Pathogen_Gene_Name", "Matched_seq", "Match_Length",
        "Pathogen_Start", "Pathogen_End",
        "Epitope_Start", "Epitope_End"
    ])

    return match_df, total_matches


match_df, total_matches = find_matches(pathogen_data, A)

print(f"Total raw matches: {total_matches}")
print(f"Before overlap filtering: {len(match_df)}")

# keep only non-overlapping matches within each protein
def keep_non_overlapping_hits(group: pd.DataFrame) -> pd.DataFrame:
    """
    Keeps longest matches and removes overlapping shorter ones
    within each protein.
    """

    group = group.sort_values(
        ["Match_Length", "Pathogen_Start", "Pathogen_End"],
        ascending=[False, True, True]
    )

    kept_rows = []
    kept_intervals = []

    for _, row in group.iterrows():
        start = row["Pathogen_Start"]
        end = row["Pathogen_End"]

        overlaps = any(
            not (end < ks or start > ke)
            for ks, ke in kept_intervals
        )

        if not overlaps:
            kept_rows.append(row)
            kept_intervals.append((start, end))

    return pd.DataFrame(kept_rows)


match_df = (
    match_df
    .groupby(
        ["Assay_ID", "Pathogen_Protein_ID", "Strain"],
        dropna=False,
        group_keys=False
    )
    .apply(keep_non_overlapping_hits)
    .reset_index(drop=True)
)

print(f"After non-overlapping filtering: {len(match_df)}")

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

full_result = all_epitopes.merge(
    match_df,
    how="left",
    on=["Assay_ID", "Epitope_Source", "Disease", "IEDB_Protein_ID"]
)

# Fill missing values
full_result["Match_Length"] = full_result["Match_Length"].fillna(
    full_result["Sequence"].str.len()
)

full_result["Pathogen_End"] = (
    full_result["Pathogen_Start"] + full_result["Match_Length"] - 1
)

full_result["Matched"] = ~full_result["Pathogen_Protein_ID"].isna()

# save results
full_result.to_csv(output_path, index=False)

print(f"Saved results to: {output_path}")
# %%
