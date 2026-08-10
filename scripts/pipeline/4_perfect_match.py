from collections import defaultdict

import ahocorasick
import pandas as pd
from tqdm import tqdm

from crossreactivity.config import config
from crossreactivity.io import DATA_DIR, load_csv, save_csv


def load_input_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    print("Loading data...")

    protein_sequences = load_csv(
        DATA_DIR / "intermediate/protein_sequences.csv"
    )
    protein_occurrences = load_csv(
        DATA_DIR / "intermediate/protein_occurrences.csv"
    )
    iedb_data = load_csv(DATA_DIR / "intermediate/wrangled_IEDB.csv")

    sequence_cols = {"Protein_ID", "Sequence", "Annotation", "Gene_name"}
    occurrence_cols = {
        "Protein_ID",
        "Proteome_ID",
        "Scientific_name",
        "Genus_species",
        "Strain",
    }

    for name, data, required in [
        ("protein_sequences.csv", protein_sequences, sequence_cols),
        ("protein_occurrences.csv", protein_occurrences, occurrence_cols),
    ]:
        missing = required.difference(data.columns)
        if missing:
            raise ValueError(f"{name} is missing columns: {sorted(missing)}")

    if protein_sequences["Protein_ID"].duplicated().any():
        raise ValueError(
            "protein_sequences.csv must contain exactly one row per Protein_ID."
        )

    if protein_occurrences.duplicated(["Protein_ID", "Proteome_ID"]).any():
        raise ValueError(
            "protein_occurrences.csv contains duplicate "
            "(Protein_ID, Proteome_ID) pairs."
        )

    for name, data, columns in [
        ("protein_sequences.csv", protein_sequences, ["Protein_ID", "Sequence"]),
        (
            "protein_occurrences.csv",
            protein_occurrences,
            ["Protein_ID", "Proteome_ID", "Genus_species"],
        ),
    ]:
        for column in columns:
            missing = data[column].isna() | data[column].astype(str).str.strip().eq("")
            if missing.any():
                raise ValueError(
                    f"{name} contains {missing.sum():,} missing {column} values."
                )

    unknown_ids = protein_occurrences.loc[
        ~protein_occurrences["Protein_ID"].isin(protein_sequences["Protein_ID"]),
        "Protein_ID",
    ].unique()
    if len(unknown_ids):
        raise ValueError(
            "Occurrence Protein_IDs missing from protein_sequences.csv. "
            f"Examples: {unknown_ids[:10].tolist()}"
        )

    return protein_sequences, protein_occurrences, iedb_data


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
    protein_sequences: pd.DataFrame,
    automaton: ahocorasick.Automaton,
) -> pd.DataFrame:
    """
    Find all exact substring matches in pathogen protein sequences.

    Each unique Protein_ID is searched once. Occurrence metadata is attached
    only after duplicate, nested-match, and IEDB-region filtering.
    """

    matches = []

    print("Matching sequences...")

    for row in tqdm(
        protein_sequences.itertuples(index=False),
        total=len(protein_sequences),
    ):
        seq = row.Sequence
        pid = row.Protein_ID
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

    Occurrence metadata is added afterward, so the final long table preserves
    disease, assay, epitope, pathogen, match, and region traceability.
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


def expand_pathogen_metadata(
    match_df: pd.DataFrame,
    protein_occurrences: pd.DataFrame,
) -> pd.DataFrame:
    """
    Expand retained protein-level matches to one row per proteome occurrence.
    """
    if match_df.empty:
        return match_df.copy()

    print("Expanding retained matches to organism-level rows...")

    metadata_columns = [
        "Proteome_ID",
        "Pathogen_Organism",
        "Pathogen_Scientific_name",
        "Pathogen_Strain",
    ]
    occurrences = protein_occurrences.rename(
        columns={
            "Genus_species": "Pathogen_Organism",
            "Scientific_name": "Pathogen_Scientific_name",
            "Strain": "Pathogen_Strain",
        }
    )
    original_columns = match_df.columns.tolist()
    expanded_df = match_df.merge(
        occurrences[["Protein_ID", *metadata_columns]],
        left_on="Pathogen_Protein_ID",
        right_on="Protein_ID",
        how="left",
        validate="many_to_many",
        indicator=True,
    )

    unmatched = expanded_df["_merge"].eq("left_only")
    if unmatched.any():
        examples = expanded_df.loc[
            unmatched, "Pathogen_Protein_ID"
        ].drop_duplicates().head(10).tolist()
        raise ValueError(
            "Matched proteins without occurrence metadata. "
            f"Examples: {examples}"
        )

    expanded_df = expanded_df.drop(columns=["Protein_ID", "_merge"])
    insert_at = original_columns.index("Pathogen_Protein_ID") + 1
    output_columns = (
        original_columns[:insert_at]
        + metadata_columns
        + original_columns[insert_at:]
    )
    expanded_df = expanded_df[output_columns]

    print(f"Retained protein-level matches: {len(match_df)}")
    print(f"Organism-level match rows after expansion: {len(expanded_df)}")

    return expanded_df


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
    protein_sequences, protein_occurrences, iedb_data = load_input_data()

    epitopes = prepare_epitopes(iedb_data)

    print("Building automaton...")
    automaton = build_automaton(epitopes)

    match_df = find_matches(
        protein_sequences=protein_sequences,
        automaton=automaton,
    )

    print(f"Raw matches from Aho-Corasick: {len(match_df)}")

    match_df = remove_exact_duplicate_matches(match_df)

    match_df = remove_nested(match_df)

    match_df_with_regions = add_iedb_region_ids(match_df)

    match_df_with_regions = expand_pathogen_metadata(
        match_df_with_regions,
        protein_occurrences,
    )

    long_output_path = DATA_DIR / "proccesed/iedb_match_regions_long.csv"
    save_csv(match_df_with_regions, long_output_path)
    print(f"Saved region-labeled long matches to: {long_output_path}")

    region_df = summarize_regions(match_df_with_regions)

    summary_output_path = DATA_DIR / "proccesed/iedb_match_regions.csv"
    save_csv(region_df, summary_output_path)
    print(f"Saved region summary to: {summary_output_path}")


if __name__ == "__main__":
    main()