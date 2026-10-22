import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from Bio import Align, SeqIO
from Bio.Align import substitution_matrices
from tqdm import tqdm

from crossreactivity.io import load_csv, save_csv, DATA_DIR, clean_id


# ---------------------- Load data ---------------------- #

matches = load_csv(DATA_DIR / "proccesed/iedb_match_regions_long.csv")
pathogens = load_csv(DATA_DIR / "intermediate/protein_sequences.csv")
fasta_path = DATA_DIR / "intermediate/matched_iedb_source_proteins.fasta"

matches["IEDB_Protein_ID"] = matches["IEDB_Protein_ID"].map(clean_id)
matches["Pathogen_Protein_ID"] = matches["Pathogen_Protein_ID"].map(clean_id)
pathogens["Protein_ID"] = pathogens["Protein_ID"].map(clean_id)

matches = matches.dropna(
    subset=["IEDB_Protein_ID", "Pathogen_Protein_ID"]
)


# ---------------------- Load IEDB protein sequences ---------------------- #

iedb_sequences = pd.DataFrame([
    {
        "IEDB_Protein_ID": clean_id(
            record.id.split("|")[1]
            if "|" in record.id
            else record.id.split()[0]
        ),
        "IEDB_Sequence": str(record.seq),
    }
    for record in SeqIO.parse(fasta_path, "fasta")
]).drop_duplicates("IEDB_Protein_ID")

missing = set(matches["IEDB_Protein_ID"]) - set(iedb_sequences["IEDB_Protein_ID"])
if missing:
    raise ValueError(f"Missing IEDB FASTA sequences: {sorted(missing)}")


# ---------------------- Prepare unique protein pairs ---------------------- #

iedb_info = (
    matches.groupby("IEDB_Protein_ID", as_index=False)
    .agg(
        Epitope_Source=("Epitope_Source", "first"),
        N_matches=("Match_ID", "nunique"),
    )
)

protein_pairs = (
    matches[["IEDB_Protein_ID", "Pathogen_Protein_ID"]]
    .drop_duplicates()
    .merge(iedb_info, on="IEDB_Protein_ID")
    .merge(
        pathogens[["Protein_ID", "Sequence"]],
        left_on="Pathogen_Protein_ID",
        right_on="Protein_ID",
    )
    .merge(iedb_sequences, on="IEDB_Protein_ID")
)

# ---------------------- Global alignment ---------------------- #

aligner = Align.PairwiseAligner()
aligner.mode = "global"
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -5
aligner.extend_gap_score = -1


def align_proteins(row):
    alignment = aligner.align(
        row["IEDB_Sequence"],
        row["Sequence"],
    )[0]

    return pd.Series({
        "Percent_Identity":
            alignment.counts().identities / alignment.length * 100,
        "Alignment_Length": alignment.length,
    })


tqdm.pandas(desc="Aligning proteins")

metrics = protein_pairs.progress_apply(align_proteins, axis=1)
alignment_df = pd.concat([protein_pairs, metrics], axis=1)


# ---------------------- Summarise per IEDB protein ---------------------- #

summary_df = (
    alignment_df
    .groupby("IEDB_Protein_ID", as_index=False)
    .agg(
        Epitope_Source=("Epitope_Source", "first"),
        Mean_Identity=("Percent_Identity", "mean"),
        SD_Identity=("Percent_Identity", "std"),
        Mean_Alignment_Length=("Alignment_Length", "mean"),
        SD_Alignment_Length=("Alignment_Length", "std"),
        N_matches=("N_matches", "first"),
        N_protein_pairs=("Pathogen_Protein_ID", "nunique"),
    )
)

save_csv(
    summary_df,
    DATA_DIR / "proccesed/global_alignment_summary.csv",
)


# ---------------------- Plot ---------------------- #
# add line for oles formula: %id * np.sqrt(alignment_length) > 290 that shows cutoff for homology based on alignment length and identity
# maybe simply define alignment length as the shortest of the two sequences, and then plot the line for that. This will give a better idea of what is considered homologous based on the alignment length and identity.

fig, ax = plt.subplots(figsize=(9, 7))

points = ax.scatter(
    summary_df["Mean_Identity"],
    summary_df["Mean_Alignment_Length"],
    c=summary_df["N_matches"],
    norm=LogNorm(),
    cmap="viridis",
    s=60,
    alpha=0.75,
)

ax.errorbar(
    summary_df["Mean_Identity"],
    summary_df["Mean_Alignment_Length"],
    xerr=summary_df["SD_Identity"],
    yerr=summary_df["SD_Alignment_Length"],
    fmt="none",
    color="grey",
    alpha=0.3,
)

ax.axvline(30, color="red", linestyle="--", label="30% identity")

max_alignment_length = summary_df["Mean_Alignment_Length"].max()
min_identity = max(1, 290 / np.sqrt(max_alignment_length))
identity_values = np.linspace(min_identity, 100, 500)
oles_alignment_lengths = (290 / identity_values) ** 2
ax.plot(
    identity_values,
    oles_alignment_lengths,
    color="black",
    linestyle=":",
    label=r"Oles cutoff: %ID $\times \sqrt{L}$ = 290",
)

ax.set_xlabel("Mean global percent identity (%)")
ax.set_ylabel("Mean global alignment length")
ax.legend()

fig.colorbar(
    points,
    ax=ax,
    label="Number of associated matches (log10 scale)",
)
fig.tight_layout()

fig.savefig(
    DATA_DIR / "proccesed/global_alignment_summary_plot.png",
    dpi=300,
)

plt.show()