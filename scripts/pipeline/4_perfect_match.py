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
    """
    Prepare IEDB epitopes for matching.

    Epitopes without protein coordinates are removed because the downstream
    IEDB-region assignment requires epitope_start_pos and epitope_end_pos.
    """

    iedb_data = iedb_data.copy()

    iedb_data["epitope_start_pos"] = pd.to_numeric(
        iedb_data["epitope_start_pos"],
        errors="coerce"
    )

    iedb_data["epitope_end_pos"] = pd.to_numeric(
        iedb_data["epitope_end_pos"],
        errors="coerce"
    )

    before = len(iedb_data)

    iedb_data = iedb_data.dropna(
        subset=["epitope_start_pos", "epitope_end_pos"]
    ).copy()

    after = len(iedb_data)

    print(f"Removed epitopes without protein coordinates: {before - after}")

    iedb_data["epitope_start_pos"] = iedb_data["epitope_start_pos"].astype(int)
    iedb_data["epitope_end_pos"] = iedb_data["epitope_end_pos"].astype(int)

    epitopes = list(
        zip(
            iedb_data["Assay_ID"],
            iedb_data["Protein_source"],
            iedb_data["Disease"],
            iedb_data["Disease_stage"],
            iedb_data["Protein_ID"],
            iedb_data["Sequence"],
            iedb_data["epitope_start_pos"],
            iedb_data["epitope_end_pos"],
            iedb_data["Response_measured"],
            iedb_data["Effector_cell"],
        )
    )

    return epitopes


def build_automaton(epitopes: list[tuple]) -> ahocorasick.Automaton:
    """
    Build an Aho-Corasick automaton for epitope-derived substrings.

    Each epitope is split into all substrings from full epitope length down to
    config["matching"]["min_match_length"].
    """

    min_match_len = config["matching"]["min_match_length"]

    automaton = ahocorasick.Automaton()
    sub_map = defaultdict(list)

    total_kmers = 0

    for (
        assay_id,
        source,
        disease,
        disease_stage,
        prot_id,
        epitope,
        start,
        end,
        response,
        effector_cell,
    ) in epitopes:

        for length in range(len(epitope), min_match_len - 1, -1):
            for i in range(len(epitope) - length + 1):

                total_kmers += 1

                sub = epitope[i:i + length]

                sub_map[sub].append(
                    (
                        assay_id,
                        source,
                        disease,
                        disease_stage,
                        prot_id,
                        epitope,
                        sub,
                        length,
                        start,
                        end,
                        i,  # offset of substring inside the full IEDB epitope
                        response,
                        effector_cell,
                    )
                )

    for sub, values in sub_map.items():
        automaton.add_word(sub, values)

    automaton.make_automaton()

    print(f"Total k-mers added to automaton: {total_kmers}")

    return automaton


def find_matches(
    pathogen_data: pd.DataFrame,
    automaton: ahocorasick.Automaton,
) -> pd.DataFrame:
    """
    Find all exact substring matches in pathogen protein sequences.

    Output is long format:
    one row = one assay/epitope substring match to one pathogen protein.
    """

    matches = []

    print("Matching sequences...")

    for row in tqdm(
        pathogen_data.itertuples(index=False),
        total=len(pathogen_data),
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
                disease_stage,
                ep_prot_id,
                full_epitope,
                sub,
                match_len,
                ep_start,
                ep_end,
                epitope_offset,
                response,
                effector_cell,
            ) in hit_list:

                # 1-based coordinates on pathogen protein
                pathogen_match_start = end_idx - match_len + 2
                pathogen_match_end = end_idx + 1

                # 1-based coordinates on original IEDB protein
                iedb_match_start = ep_start + epitope_offset
                iedb_match_end = iedb_match_start + match_len - 1

                matches.append(
                    (
                        assay_id,
                        source,
                        disease,
                        disease_stage,
                        ep_prot_id,
                        full_epitope,
                        pid,
                        proteome_id,
                        pathogen_organism,
                        scientific_name,
                        pathogen_strain,
                        annot,
                        gene,
                        sub,
                        match_len,
                        pathogen_match_start,
                        pathogen_match_end,
                        ep_start,
                        ep_end,
                        iedb_match_start,
                        iedb_match_end,
                        response,
                        effector_cell,
                    )
                )

    match_df = pd.DataFrame(
        matches,
        columns=[
            "Assay_ID",
            "Epitope_Source",
            "Disease",
            "Disease_stage",
            "IEDB_Protein_ID",
            "IEDB_Epitope_Sequence",
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
            "IEDB_Match_Start",
            "IEDB_Match_End",
            "Response_measured",
            "Effector_cell",
        ],
    )

    if not match_df.empty:
        match_df["IEDB_Epitope_Length"] = match_df[
            "IEDB_Epitope_Sequence"
        ].str.len()

    return match_df


def remove_exact_duplicate_matches(match_df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove exact duplicate rows.
    """

    n_duplicates = match_df.duplicated().sum()

    match_df = match_df.drop_duplicates().copy()

    print(f"Exact duplicate matches removed: {n_duplicates}")
    print(f"Matches after duplicate filtering: {len(match_df)}")

    return match_df


def add_iedb_region_ids(match_df: pd.DataFrame) -> pd.DataFrame:
    """
    Add merged IEDB-region labels to the long match table.

    This preserves row-level traceability:
    each row still links disease, assay, epitope, pathogen organism,
    pathogen protein, matched sequence, and region ID.
    """

    match_df = match_df.copy()

    region_labeled_groups = []

    grouped = match_df.groupby("IEDB_Protein_ID", dropna=False, sort=False)

    print(f"Adding IEDB-region IDs for {grouped.ngroups} IEDB proteins...")

    for iedb_protein_id, group in grouped:
        group = (
            group
            .dropna(subset=["IEDB_Match_Start", "IEDB_Match_End"])
            .sort_values(["IEDB_Match_Start", "IEDB_Match_End"])
            .copy()
        )

        if group.empty:
            continue

        running_end = group["IEDB_Match_End"].cummax()
        previous_running_end = running_end.shift(fill_value=-1)

        group["Region_Number"] = (
            group["IEDB_Match_Start"] > previous_running_end + 1
        ).cumsum()

        group["IEDB_Region_Start"] = (
            group
            .groupby("Region_Number")["IEDB_Match_Start"]
            .transform("min")
        )

        group["IEDB_Region_End"] = (
            group
            .groupby("Region_Number")["IEDB_Match_End"]
            .transform("max")
        )

        group["IEDB_Region_Length"] = (
            group["IEDB_Region_End"] - group["IEDB_Region_Start"] + 1
        )

        group["IEDB_Region_ID"] = (
            group["IEDB_Protein_ID"].astype(str)
            + ":"
            + group["IEDB_Region_Start"].astype(int).astype(str)
            + "-"
            + group["IEDB_Region_End"].astype(int).astype(str)
        )

        region_labeled_groups.append(group)

    if not region_labeled_groups:
        return pd.DataFrame()

    match_df_with_regions = pd.concat(region_labeled_groups, ignore_index=True)

    print(f"Long region-labeled matches: {len(match_df_with_regions)}")

    return match_df_with_regions


def unique_list(series: pd.Series) -> list:
    """
    Return unique non-missing values while preserving order.
    """

    return list(dict.fromkeys(series.dropna().tolist()))

def remove_nested(match_df: pd.DataFrame) -> pd.DataFrame:
    """
    Remove matches that add no new positional information.

    Within each assay–IEDB protein–pathogen protein combination, a match is
    removed only when both its pathogen interval and its IEDB interval are
    fully contained within the corresponding intervals of the same retained
    longer match.

    This preserves matches that map the same pathogen region to a distinct
    position in the IEDB protein, or vice versa.
    """
    print("Matches before nested filtering:", len(match_df))

    group_cols = [
        "Assay_ID",
        "IEDB_Protein_ID",
        "Pathogen_Protein_ID",
        "Proteome_ID",
    ]

    grouped_results = []

    for _, group in match_df.groupby(
        group_cols,
        dropna=False,
        sort=False,
    ):
        group_sorted = group.sort_values(
            [
                "Match_Length",
                "Pathogen_Start",
                "Pathogen_End",
                "IEDB_Match_Start",
                "IEDB_Match_End",
            ],
            ascending=[False, True, True, True, True],
        )

        kept_rows = []
        kept_matches = []

        for _, row in group_sorted.iterrows():
            pathogen_start = row["Pathogen_Start"]
            pathogen_end = row["Pathogen_End"]
            iedb_start = row["IEDB_Match_Start"]
            iedb_end = row["IEDB_Match_End"]

            is_nested = any(
                pathogen_start >= kept["Pathogen_Start"]
                and pathogen_end <= kept["Pathogen_End"]
                and iedb_start >= kept["IEDB_Match_Start"]
                and iedb_end <= kept["IEDB_Match_End"]
                for kept in kept_matches
            )

            if is_nested:
                continue

            kept_rows.append(row)

            kept_matches.append({
                "Pathogen_Start": pathogen_start,
                "Pathogen_End": pathogen_end,
                "IEDB_Match_Start": iedb_start,
                "IEDB_Match_End": iedb_end,
            })

        grouped_results.append(pd.DataFrame(kept_rows))

    if grouped_results:
        filtered_df = pd.concat(grouped_results, ignore_index=True)
    else:
        filtered_df = match_df.iloc[0:0].copy()

    print("Matches after nested filtering:", len(filtered_df))

    return filtered_df


def summarize_regions(match_df_with_regions: pd.DataFrame) -> pd.DataFrame:
    """
    Create a region-level summary from the long region-labeled match table.

    This output is useful for overview plots, but biological traceback should
    use the long table.
    """

    if match_df_with_regions.empty:
        return pd.DataFrame()

    print("Creating IEDB-region summary...")

    list_cols = {
        "Assay_ID": "Assay_IDs",
        "Epitope_Source": "Epitope_Sources",
        "Disease": "Diseases",
        "Disease_stage": "Disease_stages",
        "Response_measured": "Response_measured",
        "Effector_cell": "Effector_cells",
        "Pathogen_Protein_ID": "Pathogen_Protein_IDs",
        "Proteome_ID": "Proteome_IDs",
        "Pathogen_Organism": "Pathogen_Organisms",
        "Pathogen_Scientific_name": "Pathogen_Scientific_names",
        "Pathogen_Strain": "Pathogen_Strains",
        "Pathogen_Annotation": "Pathogen_Annotations",
        "Pathogen_Gene_Name": "Pathogen_Gene_Names",
        "Matched_seq": "Matched_sequences",
        "Match_Length": "Match_lengths",
    }

    region_df = (
        match_df_with_regions
        .groupby(
            [
                "IEDB_Protein_ID",
                "IEDB_Region_ID",
                "IEDB_Region_Start",
                "IEDB_Region_End",
                "IEDB_Region_Length",
            ],
            dropna=False,
            sort=False,
        )
        .agg(
            Number_of_supporting_matches=("Assay_ID", "size"),
            Number_of_unique_assays=("Assay_ID", "nunique"),
            Number_of_unique_pathogen_proteins=("Pathogen_Protein_ID", "nunique"),
            Number_of_unique_proteomes=("Proteome_ID", "nunique"),
            Number_of_unique_organisms=("Pathogen_Organism", "nunique"),
            **{
                new_col: (old_col, unique_list)
                for old_col, new_col in list_cols.items()
            },
        )
        .reset_index()
    )

    print(f"Total IEDB regions after merging: {len(region_df)}")

    return region_df


def main() -> None:
    pathogen_data, iedb_data = load_input_data()

    epitopes = prepare_epitopes(iedb_data)

    print("Building automaton...")
    automaton = build_automaton(epitopes)

    match_df = find_matches(
        pathogen_data=pathogen_data,
        automaton=automaton,
    )

    print(f"Raw matches from Aho-Corasick: {len(match_df)}")

    match_df = remove_exact_duplicate_matches(match_df)

    match_df = remove_nested(match_df)

    match_df_with_regions = add_iedb_region_ids(match_df)

    long_output_path = DATA_DIR / "proccesed/iedb_match_regions_long.csv"
    save_csv(match_df_with_regions, long_output_path)
    print(f"Saved region-labeled long matches to: {long_output_path}")

    region_df = summarize_regions(match_df_with_regions)

    summary_output_path = DATA_DIR / "proccesed/iedb_match_regions.csv"
    save_csv(region_df, summary_output_path)
    print(f"Saved region summary to: {summary_output_path}")


if __name__ == "__main__":
    main()