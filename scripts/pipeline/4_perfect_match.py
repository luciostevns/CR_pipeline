from collections import defaultdict

import ahocorasick
import pandas as pd
from tqdm import tqdm

from crossreactivity.config import config
from crossreactivity.io import DATA_DIR, load_csv, save_csv


def load_input_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    print("Loading data...")

    pathogen_data = load_csv(DATA_DIR / "intermediate/protein_metadata.csv")
    iedb_data = load_csv(DATA_DIR / "intermediate/wrangled_IEDB.csv")

    return pathogen_data, iedb_data


def prepare_epitopes(iedb_data: pd.DataFrame) -> list[tuple]:
    epitopes = list(
        zip(
            iedb_data["Assay_ID"],
            iedb_data["Protein_source"],
            iedb_data["Disease"],
            iedb_data["Protein_ID"],
            iedb_data["Sequence"],
            iedb_data["epitope_start_pos"],
            iedb_data["epitope_end_pos"],
        )
    )

    return epitopes


def build_automaton(epitopes: list[tuple]) -> ahocorasick.Automaton:
    """
    Build an Aho-Corasick automaton for efficient multi-pattern matching of epitope sequences.
    Each epitope is added to the automaton with associated metadata for later retrieval.
    """

    min_match_len = config["matching"]["min_match_length"]

    automaton = ahocorasick.Automaton()
    sub_map = defaultdict(list)

    for assay_id, source, disease, prot_id, epitope, start, end in epitopes:
        for length in range(len(epitope), min_match_len - 1, -1):
            for i in range(len(epitope) - length + 1):
                sub = epitope[i:i + length]

                sub_map[sub].append(
                    (
                        assay_id,
                        source,
                        disease,
                        prot_id,
                        sub,
                        length,
                        start,
                        end,
                    )
                )

    for sub, values in sub_map.items():
        automaton.add_word(sub, values)

    automaton.make_automaton()

    return automaton


def find_matches(pathogen_data: pd.DataFrame, automaton: ahocorasick.Automaton) -> pd.DataFrame:
    """
    Use the Aho-Corasick automaton to find all matches of epitope sequences within 
    the pathogen protein sequences. Each match is recorded with associated metadata
    from both the pathogen data and the epitope data.
    """
    
    matches = []

    print("Matching sequences...")

    for row in tqdm(
        pathogen_data.itertuples(index=False),
        total=len(pathogen_data)
    ):
        seq = row.Sequence
        pid = row.Protein_ID
        proteome_id = getattr(row, "Proteome_ID", None)
        pathogen_organism = getattr(row, "Genus_species", None)
        pathogen_strain = getattr(row, "Strain", None)
        scientific_name = getattr(row, "Scientific_name", None)
        annot = row.Annotation
        gene = getattr(row, "Gene_name", None)

        for end_idx, hit_list in automaton.iter(seq):
            for (
                assay_id,
                source,
                disease,
                ep_prot_id,
                sub,
                match_len,
                ep_start,
                ep_end,
            ) in hit_list:

                match_start = end_idx - match_len + 2  # 1-based
                match_end = end_idx + 1

                matches.append(
                    (
                        assay_id,
                        source,
                        disease,
                        ep_prot_id,
                        pid,
                        proteome_id,
                        pathogen_organism,
                        scientific_name,
                        pathogen_strain,
                        annot,
                        gene,
                        sub,
                        match_len,
                        match_start,
                        match_end,
                        ep_start,
                        ep_end,
                    )
                )

    match_df = pd.DataFrame(
        matches,
        columns=[
            "Assay_ID",
            "Epitope_Source",
            "Disease",
            "IEDB_Protein_ID",
            "Pathogen_Protein_ID",
            "Proteome_ID",
            "Pathogen_Organism",
            "Pathogen_Scientific_name",
            "Pathogen_Strain",
            "Pathogen_Annotation",
            "Pathogen_Gene_Name",
            "Matched_seq",
            "Match_Length",
            "Pathogen_Start",
            "Pathogen_End",
            "Epitope_Start",
            "Epitope_End",
        ],
    )

    return match_df


def remove_exact_duplicate_matches(match_df: pd.DataFrame) -> pd.DataFrame:

    n_duplicates = match_df.duplicated().sum()

    match_df = match_df.drop_duplicates().copy()

    print(f"Exact duplicate matches removed: {n_duplicates}")
    print(f"Matches after duplicate filtering: {len(match_df)}")

    return match_df


def keep_non_overlapping_hits(group: pd.DataFrame, min_sticking_out: int) -> tuple[pd.DataFrame, int]:
    """
    Removes a match only if the shorter of the two overlapping matches
    sticks out by less than min_sticking_out amino acids.
    """
    group = group.sort_values(
        ["Match_Length", "Pathogen_Start", "Pathogen_End"],
        ascending=[False, True, True],
    )

    kept_rows = []
    kept_intervals = []
    removed_count = 0

    for _, row in group.iterrows():
        start = row["Pathogen_Start"]
        end = row["Pathogen_End"]
        length = row["Match_Length"]

        too_redundant = False

        for kept_start, kept_end, kept_len in kept_intervals:
            overlap = max(
                0,
                min(end, kept_end) - max(start, kept_start) + 1
            )

            if overlap == 0:
                continue

            shorter_len = min(length, kept_len)
            sticking_out = shorter_len - overlap

            if sticking_out < min_sticking_out:
                too_redundant = True
                break

        if not too_redundant:
            kept_rows.append(row)
            kept_intervals.append((start, end, length))
        else:
            removed_count += 1

    return pd.DataFrame(kept_rows), removed_count


def remove_overlapping_matches(match_df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove overlapping within defined groups of matches
    based on the criteria defined in keep_non_overlapping_hits.

    """

    min_sticking_out = config["matching"]["overlap_min_sticking_out"]

    print(f"Total matches before overlap filtering: {len(match_df)}")

    group_cols = [
        "Assay_ID",
        "Pathogen_Protein_ID",
        "Proteome_ID",
        "Pathogen_Organism",
    ]

    grouped_results = []
    total_removed = 0

    for keys, group in match_df.groupby(group_cols, dropna=False, sort=False):
        filtered_group, removed = keep_non_overlapping_hits(
            group=group.copy(),
            min_sticking_out=min_sticking_out,
        )

        total_removed += removed

        for col, val in zip(
            group_cols,
            keys if isinstance(keys, tuple) else (keys,)
        ):
            filtered_group[col] = val

        grouped_results.append(filtered_group)

    match_df = pd.concat(grouped_results, ignore_index=True)

    print(f"Total overlaps removed: {total_removed}")
    print(f"Total matches after overlap filtering: {len(match_df)}")

    return match_df


def add_unmatched_epitopes(match_df: pd.DataFrame, iedb_data: pd.DataFrame) -> pd.DataFrame:
    all_epitopes = (
        iedb_data[
            [
                "Assay_ID",
                "Protein_source",
                "Disease",
                "Protein_ID",
                "Sequence",
                "epitope_start_pos",
                "epitope_end_pos",
            ]
        ]
        .rename(
            columns={
                "Protein_source": "Epitope_Source",
                "Protein_ID": "IEDB_Protein_ID",
            }
        )
        .drop_duplicates()
    )

    full_result = all_epitopes.merge(
        match_df,
        how="left",
        on=[
            "Assay_ID",
            "Epitope_Source",
            "Disease",
            "IEDB_Protein_ID",
        ],
    )

    full_result["Matched"] = ~full_result["Pathogen_Protein_ID"].isna()

    full_result.loc[full_result["Matched"], "Pathogen_End"] = (
        full_result.loc[full_result["Matched"], "Pathogen_Start"]
        + full_result.loc[full_result["Matched"], "Match_Length"]
        - 1
    )

    print(full_result["Matched"].value_counts(dropna=False))

    return full_result


def main() -> None:
    pathogen_data, iedb_data = load_input_data()

    epitopes = prepare_epitopes(iedb_data)

    print("Building automaton...")
    automaton = build_automaton(epitopes)

    match_df = find_matches(
        pathogen_data=pathogen_data,
        automaton=automaton
    )

    print(f"Raw matches from Aho-Corasick: {len(match_df)}")

    match_df = remove_exact_duplicate_matches(match_df)
    match_df = remove_overlapping_matches(match_df)

    full_result = add_unmatched_epitopes(
        match_df=match_df,
        iedb_data=iedb_data
    )

    output_path = DATA_DIR / "proccesed/perfect_matches_2_0.csv"
    save_csv(full_result, output_path)

    print(f"Saved results to: {output_path}")


if __name__ == "__main__":
    main()