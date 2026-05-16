#%%
import time
import pandas as pd
from tqdm import tqdm

from helpers import (
    load_csv,
    DATA_DIR,
    clean_id,
    write_fasta_record,
    fetch_uniprot_fasta,
    parse_fasta_text,
    gram_status
)

# Data paths
match_file = DATA_DIR / "proccesed/perfect_matches_2_0.csv"
protein_meta_file = DATA_DIR / "intermediate/protein_metadata.csv"

iedb_fasta_out = DATA_DIR / "intermediate/matched_iedb_source_proteins.fasta"
unresolved_out = DATA_DIR / "intermediate/unresolved_iedb_protein_ids.csv"

# Load data
print("Loading data file...")
match_df = load_csv(match_file)
protein_meta = load_csv(protein_meta_file)

# set batch size for pathogen FASTA export
PATHOGEN_FASTA_BATCH_SIZE = 500

# Clean and extract unique matched pathogen ids
matched_pathogen_ids = (
    match_df["Pathogen_Protein_ID"]
    .map(clean_id)
    .dropna()
    .unique()
    .tolist()
)

print(f"Unique matched pathogen Protein_IDs: {len(matched_pathogen_ids)}")

# Clean Protein_IDs in protein_meta for matching
protein_meta["Protein_ID_clean"] = protein_meta["Protein_ID"].map(clean_id)

# Subset protein_meta to only matched pathogen proteins
pathogen_subset = protein_meta[
    protein_meta["Protein_ID_clean"].isin(matched_pathogen_ids)
].copy()

print(f"Matched protein-proteome entries found in protein_metadata.csv: {len(pathogen_subset)}")

pathogen_subset["Gram_status"] = pathogen_subset["Genus_species"].map(gram_status)

gram_pos_df = pathogen_subset[pathogen_subset["Gram_status"] == "positive"]
gram_neg_df = pathogen_subset[pathogen_subset["Gram_status"] == "negative"]

unknown_df = pathogen_subset[pathogen_subset["Gram_status"].isna()]

print("\nGram status counts:")
print(pathogen_subset["Gram_status"].value_counts(dropna=False))

if len(unknown_df) > 0:
    print("\nUnmapped species:")
    print(
        unknown_df["Genus_species"]
        .drop_duplicates()
        .sort_values()
        .to_string(index=False)
    )

def export_fasta_batches(df, prefix, batch_size=500):
    n_total = len(df)
    n_batches = (n_total + batch_size - 1) // batch_size

    print(f"\nExporting {prefix}: {n_total} proteins in {n_batches} batch(es)")

    for batch_idx in range(n_batches):
        start = batch_idx * batch_size
        end = min(start + batch_size, n_total)
        batch_df = df.iloc[start:end]

        batch_file = (
            DATA_DIR / "intermediate" /
            f"{prefix}_part_{batch_idx + 1:03d}.fasta"
        )

        with open(batch_file, "w", encoding="utf-8") as out:
            for row in batch_df.itertuples(index=False):
                protein_id = clean_id(row.Protein_ID)
                proteome_id = clean_id(getattr(row, "Proteome_ID", None)) or "NA"

                organism = getattr(row, "Genus_species", None)
                organism = str(organism).strip() if pd.notna(organism) and organism else "NA"

                annotation = getattr(row, "Annotation", None)
                annotation = str(annotation).strip() if pd.notna(annotation) and annotation else "NA"

                gram = getattr(row, "Gram_status", None)
                gram = str(gram).strip() if pd.notna(gram) and gram else "NA"

                sequence = getattr(row, "Sequence", None)

                if pd.isna(sequence) or not protein_id:
                    continue

                sequence = str(sequence).strip()
                if not sequence:
                    continue

                header = (
                    f"{protein_id} "
                    f"organism={organism} "
                    f"proteome={proteome_id} "
                    f"gram={gram} "
                    f"annotation={annotation}"
                )

                write_fasta_record(out, header, sequence)

        print(f"Saved {prefix} FASTA batch {batch_idx + 1}/{n_batches}: {batch_file}")

export_fasta_batches(
    gram_pos_df,
    prefix="matched_pathogen_proteins_gram_positive",
    batch_size=PATHOGEN_FASTA_BATCH_SIZE
)

export_fasta_batches(
    gram_neg_df,
    prefix="matched_pathogen_proteins_gram_negative",
    batch_size=PATHOGEN_FASTA_BATCH_SIZE
)

print("Preparing matched IEDB source proteins...")

# Clean and extract unique matched IEDB source protein ids
matched_iedb_ids = (
    match_df["IEDB_Protein_ID"]
    .map(clean_id)
    .dropna()
    .unique()
    .tolist()
)

print(f"Unique matched IEDB source protein IDs: {len(matched_iedb_ids)}")

# vars for tracking not found proteins
resolved = 0
unresolved = []

# Fetch FASTA sequences for matched IEDB source proteins and write to output FASTA
with open(iedb_fasta_out, "w", encoding="utf-8") as out:
    for protein_id in tqdm(matched_iedb_ids, desc="Fetching IEDB proteins", unit="protein"):
        fasta_text = fetch_uniprot_fasta(protein_id)
        parsed = parse_fasta_text(fasta_text) if fasta_text else None

        if parsed is None:
            unresolved.append({"IEDB_Protein_ID": protein_id})
            continue

        header, sequence = parsed
        write_fasta_record(out, header, sequence)
        resolved += 1

        time.sleep(0.2)

print(f"Saved IEDB source protein FASTA to: {iedb_fasta_out}")
print(f"Resolved IEDB source proteins: {resolved}")
print(f"Unresolved IEDB source proteins: {len(unresolved)}")

if unresolved:
    pd.DataFrame(unresolved).to_csv(unresolved_out, index=False)
    print(f"Saved unresolved IDs to: {unresolved_out}")

print("Done.")
#%%