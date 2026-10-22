import pandas as pd
import matplotlib.pyplot as plt
from crossreactivity.io import DATA_DIR

path = DATA_DIR / "raw/1229370_NetMHCIIpan.xls"
with open(path) as f:
    header = next(i for i, line in enumerate(f) if "Peptide" in line)

pred = pd.read_csv(path, sep="\t", header=header)
pred.columns = pred.columns.str.lower().str.replace(r"[^a-z0-9]", "", regex=True)

regions = pd.read_csv(
    DATA_DIR / "proccesed/iedb_match_regions_long_metadata_annotated.csv"
)

# Match exact epitope sequences, then retain cores present in Matched_seq.
pairs = regions.reset_index(names="region_row").merge(
    pred, left_on="IEDB_Epitope_Sequence", right_on="peptide"
)
pairs = pairs[
    pairs.apply(
        lambda r: isinstance(r["core"], str)
        and bool(r["core"])
        and isinstance(r["Matched_seq"], str)
        and r["core"] in r["Matched_seq"],
        axis=1,
    )
]

matched = regions.loc[pairs["region_row"].unique()].copy()
scores = pairs.drop_duplicates(["peptide", "mhc", "core"]).copy()
scores[["rankba", "rankel"]] = scores[["rankba", "rankel"]].apply(pd.to_numeric)

print(f"Retained rows: {len(matched):,}")
print(f"Unique Match_IDs: {matched['Match_ID'].nunique():,}")
print(f"Unique peptide–MHC–core predictions: {len(scores):,}")
print(scores[["rankba", "rankel"]].agg(["count", "nunique", "min", "max"]))

print(
    scores[["peptide", "mhc", "core", "rankba", "rankel"]]
    .to_string(index=False)
)

print(
    pairs.groupby(["peptide", "mhc", "core"])["Match_ID"]
    .nunique()
    .rename("unique_matches")
    .to_string()
)

print(matched["Assay_ID"].dropna().drop_duplicates().to_string(index=False))