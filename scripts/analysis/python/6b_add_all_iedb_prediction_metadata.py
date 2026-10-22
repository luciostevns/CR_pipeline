"""Add all-IEDB DeepLoc and NetSurfP predictions to wrangled_IEDB.csv.

Inputs:
  data/raw/netsurfp/all_iedb_protein_netsurfp_predictions.csv
  data/raw/deeploc/all_iedb_protein_deeploc_predictions.csv
Outputs:
  data/proccesed/iedb_data_with_metadata.csv
  data/proccesed/iedb_prediction_coverage.csv

Preserves every original IEDB row. RSA uses inclusive epitope coordinates,
not the shorter matched region. Partial residue coverage is explicitly flagged.
"""

import re
from bisect import bisect_left, bisect_right
from functools import lru_cache
from itertools import accumulate

import pandas as pd

from crossreactivity.io import DATA_DIR, load_csv, save_csv


def normalize_id(value):
    if pd.isna(value) or not str(value).strip():
        return None
    token = str(value).strip().lstrip(">").split()[0].rsplit("/", 1)[-1]
    if "|" in token:
        token = token.split("|")[1]
    else:
        legacy = re.match(r"(?:sp|tr)_([^_]+)_", token, flags=re.I)
        if legacy:
            token = legacy.group(1)
    return token.strip().upper() or None


def load_locations(path):
    data = pd.read_csv(path, dtype=str)
    data.columns = data.columns.str.strip()
    required = {"Protein_ID", "Localizations"}
    if not required.issubset(data.columns):
        raise ValueError(f"DeepLoc CSV requires columns: {sorted(required)}")
    data["key"] = data["Protein_ID"].map(normalize_id)
    data["location"] = data["Localizations"].astype("string").str.strip().replace("", pd.NA)
    data = data.dropna(subset=["key", "location"])[["key", "location"]].drop_duplicates()
    conflicts = data.groupby("key")["location"].nunique()
    if conflicts.gt(1).any():
        raise ValueError(f"Conflicting DeepLoc locations: {conflicts[conflicts.gt(1)].index.tolist()}")
    return data.set_index("key")["location"].to_dict()


def load_rsa(path):
    # Same 21-column NetSurfP format as step 6. Splitting from the right
    # also supports older FASTA headers containing commas.
    rows = []
    with path.open(encoding="utf-8-sig") as handle:
        next(handle, None)
        for line_number, line in enumerate(handle, 2):
            if not line.strip():
                continue
            fields = line.strip().rsplit(",", 20)
            if len(fields) != 21:
                raise ValueError(f"Malformed NetSurfP row: {path}:{line_number}")
            rows.append((normalize_id(fields[0].strip('"')), fields[2], fields[3]))
    data = pd.DataFrame(rows, columns=["key", "position", "rsa"])
    for column in ["position", "rsa"]:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data = data.dropna().copy()
    valid = (data["position"].ge(1) & data["position"].mod(1).eq(0)
             & data["rsa"].between(0, 1))
    if not valid.all():
        raise ValueError("NetSurfP contains invalid residue positions or RSA outside 0–1.")
    data["position"] = data["position"].astype(int)
    data = data.drop_duplicates()
    if data.duplicated(["key", "position"]).any():
        raise ValueError("Conflicting NetSurfP RSA values for the same protein/residue.")
    lookup = {}
    for protein_id, group in data.groupby("key"):
        group = group.sort_values("position")
        lookup[protein_id] = (
            group["position"].tolist(),
            [0.0, *accumulate(group["rsa"].tolist())],
        )
    return lookup


def make_region_summarizer(lookup):
    @lru_cache(maxsize=None)
    def summarize(protein_id, start, end):
        if protein_id is None:
            return None, 0, None, "missing_protein_id"
        if start is None or end is None:
            return None, 0, None, "missing_coordinates"
        if not (start >= 1 and end >= 1 and float(start).is_integer() and float(end).is_integer()):
            return None, 0, None, "invalid_coordinates"
        start, end = sorted((int(start), int(end)))
        expected = end - start + 1
        if protein_id not in lookup:
            return None, 0, 0.0, "missing_prediction"
        positions, sums = lookup[protein_id]
        left = bisect_left(positions, start)
        right = bisect_right(positions, end)
        count = right - left
        if not count:
            return None, 0, 0.0, "no_residues_in_interval"
        # Like step 6, average the available residues; flag partial coverage.
        return ((sums[right] - sums[left]) / count, count, count / expected,
                "complete" if count == expected else "partial")
    return summarize


def add_metadata(iedb_data, locations, rsa_lookup):
    result = iedb_data.copy()
    result["IEDB_Protein_ID_clean"] = result["Protein_ID"].map(normalize_id)
    result["IEDB_prot_location"] = result["IEDB_Protein_ID_clean"].map(locations)
    starts = pd.to_numeric(result["epitope_start_pos"], errors="coerce")
    ends = pd.to_numeric(result["epitope_end_pos"], errors="coerce")
    summarize = make_region_summarizer(rsa_lookup)
    values = [summarize(pid, None if pd.isna(start) else start,
                        None if pd.isna(end) else end)
              for pid, start, end in zip(result["IEDB_Protein_ID_clean"], starts, ends)]
    columns = ["rsa_epitope_mean", "rsa_epitope_n_residues",
               "rsa_epitope_coverage", "rsa_epitope_status"]
    for offset, column in enumerate(columns):
        result[column] = [value[offset] for value in values]
    return result


def main():
    iedb_data = load_csv(DATA_DIR / "intermediate/wrangled_IEDB.csv")
    locations = load_locations(
        DATA_DIR / "raw/deeploc/all_iedb_protein_deeploc_predictions.csv"
    )
    rsa_lookup = load_rsa(
        DATA_DIR / "raw/netsurfp/all_iedb_protein_netsurfp_predictions.csv"
    )
    result = add_metadata(iedb_data, locations, rsa_lookup)
    output = DATA_DIR / "proccesed/iedb_data_with_metadata.csv"
    save_csv(result, output)

    coverage = result[["IEDB_Protein_ID_clean"]].dropna().drop_duplicates()
    coverage["has_location"] = coverage["IEDB_Protein_ID_clean"].isin(locations)
    coverage["has_rsa_prediction"] = coverage["IEDB_Protein_ID_clean"].isin(rsa_lookup)
    save_csv(coverage, DATA_DIR / "proccesed/iedb_prediction_coverage.csv")
    print(f"Saved all {len(result):,} IEDB rows to {output}")
    print(f"Rows with location: {result['IEDB_prot_location'].notna().sum():,}")
    print("Epitope RSA coverage:")
    print(result["rsa_epitope_status"].value_counts().to_string())


if __name__ == "__main__":
    main()
