"""Optional IEDB coverage analysis; run with --plot to also save a figure.

Uses step 7's annotated matches, wrangled IEDB data, and step 5's FASTA.
Coordinates are 1-based and inclusive. Output groups match the QMD:
IEDB_Protein_ID, Epitope_Source, Disease. Only matched groups are reported.
"""

import argparse
import json
import re

import pandas as pd
from Bio import SeqIO

from crossreactivity.io import DATA_DIR, load_csv, save_csv


GROUP_COLS = ["IEDB_Protein_ID", "Epitope_Source", "Disease"]
HISTONE_CONTAMINANT_IDS = {"A0AB35K8W1", "A0AB35K7K7", "A0AB35K6S2"}
CATEGORIES = {
    "matched_pct": "Matched",
    "studied_but_not_matched_pct": "Studied but not matched",
    "not_studied_pct": "Not studied",
}


def parse_uniprot_id(value):
    if pd.isna(value):
        return None
    text = str(value).strip()
    if not text:
        return ""
    token = re.sub(r"^.*/", "", text.split()[0])
    return (token.split("|")[1] if "|" in token else token).strip().upper()


def merge_intervals(intervals):
    """Merge overlapping or directly adjacent inclusive intervals."""
    normalized = sorted({
        (int(min(start, end)), int(max(start, end)))
        for start, end in intervals
        if pd.notna(start) and pd.notna(end)
    })
    merged = []
    for start, end in normalized:
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def intersect_intervals(a, b):
    """Intersect sorted, merged intervals without an all-pairs expansion."""
    intersections = []
    i = j = 0
    while i < len(a) and j < len(b):
        start = max(a[i][0], b[j][0])
        end = min(a[i][1], b[j][1])
        if start <= end:
            intersections.append((start, end))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return merge_intervals(intersections)


def interval_length(intervals):
    return sum(end - start + 1 for start, end in intervals)


def grouped_intervals(data, start_col, end_col):
    result = {}
    for key, group in data.groupby(GROUP_COLS, dropna=False, sort=True):
        # None provides a stable dictionary key for missing source/disease.
        key = tuple(None if pd.isna(value) else value for value in key)
        result[key] = merge_intervals(zip(group[start_col], group[end_col]))
    return result


def calculate_coverage(matches, iedb_data, protein_lengths):
    matches = matches.copy()
    matches["IEDB_Protein_ID"] = matches["IEDB_Protein_ID"].map(parse_uniprot_id)
    matches["Pathogen_Protein_ID"] = (
        matches["Pathogen_Protein_ID"].astype("string").str.strip().str.upper()
    )
    contaminant = (
        matches["Pathogen_Protein_ID"].isin(HISTONE_CONTAMINANT_IDS)
        & matches["Epitope_Source"].fillna("").str.contains("histone", case=False)
    )
    matches = matches.loc[~contaminant].copy()
    matches = matches[GROUP_COLS + [
        "Match_ID", "IEDB_Match_Start", "IEDB_Match_End"
    ]].drop_duplicates()

    studied = iedb_data.rename(columns={
        "Protein_ID": "IEDB_Protein_ID", "Protein_source": "Epitope_Source"
    }).copy()
    studied["IEDB_Protein_ID"] = studied["IEDB_Protein_ID"].map(parse_uniprot_id)
    for data, columns in [
        (studied, ["epitope_start_pos", "epitope_end_pos"]),
        (matches, ["IEDB_Match_Start", "IEDB_Match_End"]),
    ]:
        for column in columns:
            data[column] = pd.to_numeric(data[column], errors="raise")
    studied = studied.dropna(subset=[
        "IEDB_Protein_ID", "epitope_start_pos", "epitope_end_pos"
    ])
    studied_groups = grouped_intervals(studied, "epitope_start_pos", "epitope_end_pos")
    matched_groups = grouped_intervals(matches, "IEDB_Match_Start", "IEDB_Match_End")

    rows = []
    for key, matched in matched_groups.items():
        if key not in studied_groups:
            raise ValueError(f"Matched group cannot be linked to IEDB data: {key}")
        protein_id, source, disease = key
        length = protein_lengths.get(protein_id)
        if length is None or length <= 0:
            raise ValueError(f"Missing or invalid FASTA length: {protein_id}")
        studied_intervals = studied_groups[key]
        if any(start < 1 or end > length for start, end in studied_intervals + matched):
            raise ValueError(f"Coordinates outside protein length {length}: {key}")

        intersection = intersect_intervals(studied_intervals, matched)
        studied_length = interval_length(studied_intervals)
        matched_length = interval_length(intersection)
        unmatched_length = studied_length - matched_length
        not_studied_length = length - studied_length
        display = ["NA" if value is None else str(value) for value in key]
        rows.append({
            **dict(zip(GROUP_COLS, key)),
            "Protein_Length": length,
            "studied_length": studied_length,
            "matched_within_studied_length": matched_length,
            "studied_but_not_matched_length": unmatched_length,
            "not_studied_length": not_studied_length,
            "matched_pct": 100 * matched_length / length,
            "studied_but_not_matched_pct": 100 * unmatched_length / length,
            "not_studied_pct": 100 * not_studied_length / length,
            "Protein_Plot_ID": "___".join(display),
            "Protein_Label": f"{display[1]}\n({display[2]})",
            "studied_intervals": json.dumps(studied_intervals),
            "matched_intervals": json.dumps(matched),
            "matched_within_studied_intervals": json.dumps(intersection),
        })
    columns = GROUP_COLS + [
        "Protein_Length", "studied_length", "matched_within_studied_length",
        "studied_but_not_matched_length", "not_studied_length", *CATEGORIES,
        "Protein_Plot_ID", "Protein_Label", "studied_intervals",
        "matched_intervals", "matched_within_studied_intervals",
    ]
    return pd.DataFrame(rows, columns=columns)


def plot_coverage(coverage, output_path):
    import matplotlib.pyplot as plt

    if coverage.empty:
        print("No matched groups to plot.")
        return
    ordered = coverage.sort_values(["Protein_Length", "Protein_Plot_ID"])
    fig, ax = plt.subplots(figsize=(12, max(5, len(ordered) * 0.23)))
    left = pd.Series(0.0, index=ordered.index)
    colors = {"Matched": "#F8766D", "Studied but not matched": "#00BA38", "Not studied": "#619CFF"}
    # ggplot's default stack puts the last factor level nearest zero.
    for column, label in reversed(list(CATEGORIES.items())):
        ax.barh(range(len(ordered)), ordered[column], left=left,
                height=0.8, label=label, color=colors[label])
        left = left + ordered[column]
    ax.set_yticks(range(len(ordered)), ordered["Protein_Label"], fontsize=5)
    ax.set_xlim(0, 100)
    ax.set_xlabel("Percentage of full IEDB source-protein length")
    ax.set_ylabel("IEDB source protein")
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], title="Region type",
              loc="upper left", bbox_to_anchor=(1, 1))
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plot", action="store_true", help="Also save the coverage plot")
    args = parser.parse_args()
    processed = DATA_DIR / "proccesed"
    intermediate = DATA_DIR / "intermediate"
    matches = load_csv(processed / "iedb_match_regions_long_metadata_annotated.csv")
    iedb_data = load_csv(intermediate / "wrangled_IEDB.csv")
    lengths = {}
    for record in SeqIO.parse(intermediate / "matched_iedb_source_proteins.fasta", "fasta"):
        protein_id = parse_uniprot_id(record.id)
        length = len(record.seq)
        if protein_id in lengths and lengths[protein_id] != length:
            raise ValueError(f"Conflicting FASTA lengths: {protein_id}")
        lengths[protein_id] = length

    coverage = calculate_coverage(matches, iedb_data, lengths)
    long = coverage.melt(
        id_vars=GROUP_COLS + ["Protein_Plot_ID", "Protein_Label", "Protein_Length"],
        value_vars=list(CATEGORIES), var_name="region_type", value_name="percentage",
    )
    long["region_type"] = long["region_type"].map(CATEGORIES)
    save_csv(coverage, processed / "iedb_protein_coverage.csv")
    save_csv(long, processed / "iedb_protein_coverage_long.csv")
    if args.plot:
        plot_coverage(coverage, processed / "iedb_protein_coverage.png")
    print(f"Saved coverage for {len(coverage):,} protein/source/disease groups.")


if __name__ == "__main__":
    main()
