import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from Bio import SeqIO
from matplotlib.lines import Line2D
from pycodamath.extra import norm

import pycodamath as coda

from crossreactivity.io import DATA_DIR, load_csv, clean_id


match_df = load_csv(DATA_DIR / "proccesed/iedb_match_regions_long.csv")
match_df["IEDB_Protein_ID"] = match_df["IEDB_Protein_ID"].map(clean_id)

# Remove proteins whose annotation contains "histone", regardless of capitalization
match_df = match_df[
    ~match_df["Epitope_Source"]
    .fillna("")
    .str.contains("histone", case=False, regex=False)
].copy()


lengths = {
    clean_id(
        record.id.split("|")[1]
        if "|" in record.id
        else record.id.split()[0]
    ): len(record.seq)
    for record in SeqIO.parse(
        DATA_DIR / "intermediate/matched_iedb_source_proteins.fasta",
        "fasta",
    )
}


def count_regions(group):
    studied = set().union(*[
        range(int(start), int(end) + 1)
        for start, end in zip(
            group["Epitope_Start"],
            group["Epitope_End"],
        )
        if pd.notna(start) and pd.notna(end)
    ])

    matched = set().union(*[
        range(int(start), int(end) + 1)
        for start, end in zip(
            group["IEDB_Match_Start"],
            group["IEDB_Match_End"],
        )
        if pd.notna(start) and pd.notna(end)
    ]) & studied

    return pd.Series({
        "Matched": len(matched),
        "Studied but not matched": len(studied - matched),
        "Not studied": lengths[group.name] - len(studied),
    })


count_matrix = (
    match_df
    .groupby("IEDB_Protein_ID")
    .apply(count_regions, include_groups=False)
    .reset_index()
)


percentage_matrix = (
    count_matrix
    .set_index("IEDB_Protein_ID")
    .coda.closure(100)
)

