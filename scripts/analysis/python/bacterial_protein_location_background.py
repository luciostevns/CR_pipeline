"""Export a reproducible random background sample for DeepLoc prediction."""

import pandas as pd

from crossreactivity.io import DATA_DIR, load_csv, clean_id, save_csv
from crossreactivity.fasta import write_fasta_record
from crossreactivity.reference_data import GRAM_STATUS


SAMPLE_SIZE = 5000
RANDOM_SEED = 42
BATCH_SIZE = 500
MIN_LENGTH = 20
MAX_LENGTH = 6000
OUTPUT_DIR = DATA_DIR / "intermediate/background_prediction"


def prepare_background():
    proteins = load_csv(DATA_DIR / "intermediate/protein_sequences.csv")
    proteins["Protein_ID"] = proteins["Protein_ID"].map(clean_id)
    if proteins["Protein_ID"].isna().any():
        raise ValueError("Missing Protein_ID in protein_sequences.csv")
    if proteins["Protein_ID"].duplicated().any():
        raise ValueError("Expected one sequence row per Protein_ID")
    # Remove sequence whitespace before checking the predictor's length limits.
    proteins["Sequence"] = proteins["Sequence"].astype("string").str.replace(r"\s+", "", regex=True)
    lengths = proteins["Sequence"].str.len()
    eligible = lengths.between(MIN_LENGTH, MAX_LENGTH).fillna(False)
    print(f"Proteins within {MIN_LENGTH}–{MAX_LENGTH} amino acids: {eligible.sum():,} / {len(proteins):,}")
    proteins = proteins.loc[eligible].copy()
    if SAMPLE_SIZE > len(proteins):
        raise ValueError("Sample size exceeds the number of proteins within the length limits")

    # Sample before joining occurrences: every unique protein has equal chance.
    sample = proteins.sample(n=SAMPLE_SIZE, random_state=RANDOM_SEED).copy()
    if sample["Sequence"].isna().any() or sample["Sequence"].str.strip().eq("").any():
        raise ValueError("Sample contains missing sequences")

    occurrences = load_csv(DATA_DIR / "intermediate/protein_occurrences.csv")
    occurrences["Protein_ID"] = occurrences["Protein_ID"].map(clean_id)
    occurrences = occurrences[occurrences["Protein_ID"].isin(sample["Protein_ID"])].copy()
    occurrences["Genus_species"] = occurrences["Genus_species"].astype("string").str.strip()
    occurrences["Gram_status"] = occurrences["Genus_species"].map(GRAM_STATUS)

    missing = occurrences.loc[occurrences["Gram_status"].isna(), "Genus_species"]
    if not missing.empty:
        raise ValueError(
            "Missing GRAM_STATUS mapping for: "
            + ", ".join(missing.fillna("<missing species>").unique())
        )
    if not set(occurrences["Gram_status"]).issubset({"positive", "negative"}):
        raise ValueError("GRAM_STATUS must contain 'positive' or 'negative'")
    absent = set(sample["Protein_ID"]) - set(occurrences["Protein_ID"])
    if absent:
        raise ValueError(f"Sampled proteins without occurrences: {sorted(absent)[:10]}")

    # Match the original exporter: one input per protein/Gram class, not strain.
    membership = occurrences[["Protein_ID", "Gram_status"]].drop_duplicates()
    predictions = membership.merge(sample, on="Protein_ID", validate="many_to_one")
    suffix = predictions["Gram_status"].map({"positive": "GP", "negative": "GN"})
    predictions["Prediction_ID"] = predictions["Protein_ID"] + "_" + suffix
    return sample, occurrences, predictions


def export_batches(predictions):
    manifests = []
    for gram in ("positive", "negative"):
        group = predictions[predictions["Gram_status"].eq(gram)]
        for start in range(0, len(group), BATCH_SIZE):
            batch = group.iloc[start:start + BATCH_SIZE].copy()
            filename = (
                f"background_pathogen_proteins_deeploc_gram_{gram}"
                f"_part_{start // BATCH_SIZE + 1:03d}.fasta"
            )
            with open(OUTPUT_DIR / filename, "w", encoding="utf-8") as out:
                for row in batch.itertuples(index=False):
                    annotation = getattr(row, "Annotation", "NA")
                    annotation = "NA" if pd.isna(annotation) else " ".join(str(annotation).split())
                    write_fasta_record(
                        out,
                        f"{row.Prediction_ID} annotation={annotation}",
                        str(row.Sequence).strip(),
                    )
            batch["FASTA_file"] = filename
            manifests.append(batch[["Prediction_ID", "Protein_ID", "Gram_status", "FASTA_file"]])
            print(f"Saved {filename}: {len(batch)} proteins")
    return pd.concat(manifests, ignore_index=True)


def main():
    sample, occurrences, predictions = prepare_background()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    # Prevent old batches being mistaken for part of a new sample.
    if any(OUTPUT_DIR.glob("*.fasta")):
        raise FileExistsError(f"FASTA files already exist in {OUTPUT_DIR}; use a new output directory")
    save_csv(sample, OUTPUT_DIR / "background_sample.csv")
    save_csv(occurrences, OUTPUT_DIR / "background_sample_occurrences.csv")
    manifest = export_batches(predictions)
    save_csv(manifest, OUTPUT_DIR / "background_deeploc_prediction_manifest.csv")
    print(f"\nRandom seed: {RANDOM_SEED}")
    print(f"Unique sampled proteins: {len(sample):,}")
    print(f"Prediction inputs: {len(predictions):,}")
    print(predictions["Gram_status"].value_counts().to_string())
    print(f"Output directory: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
