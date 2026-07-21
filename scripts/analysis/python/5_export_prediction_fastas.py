import time

import pandas as pd
from tqdm import tqdm

from crossreactivity.io import DATA_DIR, load_csv, clean_id, save_csv
from crossreactivity.fasta import write_fasta_record, parse_fasta_text
from crossreactivity.uniprot import fetch_uniprot_fasta
from crossreactivity.reference_data import GRAM_STATUS


PATHOGEN_FASTA_BATCH_SIZE = 500


def load_input_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    match_file = DATA_DIR / "proccesed/iedb_match_regions_long.csv"
    protein_meta_file = DATA_DIR / "intermediate/protein_metadata.csv"

    print("Loading data files...")
    match_df = load_csv(match_file)
    protein_meta = load_csv(protein_meta_file)

    return match_df, protein_meta


def get_matched_pathogen_subset(
    match_df: pd.DataFrame,
    protein_meta: pd.DataFrame
) -> pd.DataFrame:
    matched_pathogen_ids = (
        match_df["Pathogen_Protein_ID"]
        .map(clean_id)
        .dropna()
        .unique()
        .tolist()
    )

    print(f"Unique matched pathogen Protein_IDs: {len(matched_pathogen_ids)}")

    protein_meta = protein_meta.copy()
    protein_meta["Protein_ID_clean"] = protein_meta["Protein_ID"].map(clean_id)

    pathogen_subset = protein_meta[
        protein_meta["Protein_ID_clean"].isin(matched_pathogen_ids)
    ].copy()

    print(
        "Matched protein-proteome entries found in protein_metadata.csv: "
        f"{len(pathogen_subset)}"
    )

    sequence_check = (
        pathogen_subset
        .dropna(subset=["Protein_ID_clean", "Sequence"])
        .groupby("Protein_ID_clean")["Sequence"]
        .nunique()
    )

    conflicting_sequences = sequence_check[sequence_check > 1]

    print(
        "Protein IDs with more than one unique sequence: "
        f"{len(conflicting_sequences)}"
    )

    if len(conflicting_sequences) > 0:
        print("Examples of conflicting protein IDs:")
        print(conflicting_sequences.head(20))

    pathogen_subset = (
        pathogen_subset
        .sort_values(["Protein_ID_clean", "Proteome_ID"])
        .drop_duplicates(subset=["Protein_ID_clean"], keep="first")
        .copy()
    )

    print(
        "Unique matched pathogen proteins retained for FASTA export: "
        f"{len(pathogen_subset)}"
    )

    return pathogen_subset


def add_gram_status(pathogen_subset: pd.DataFrame) -> pd.DataFrame:
    pathogen_subset = pathogen_subset.copy()
    pathogen_subset["Gram_status"] = pathogen_subset["Genus_species"].map(GRAM_STATUS)

    print("\nGram status counts:")
    print(pathogen_subset["Gram_status"].value_counts(dropna=False))

    unknown_df = pathogen_subset[pathogen_subset["Gram_status"].isna()]

    if len(unknown_df) > 0:
        print("\nUnmapped species:")
        print(
            unknown_df["Genus_species"]
            .drop_duplicates()
            .sort_values()
            .to_string(index=False)
        )

    return pathogen_subset


def export_fasta_batches(
    df: pd.DataFrame,
    prefix: str,
    batch_size: int = PATHOGEN_FASTA_BATCH_SIZE
) -> None:
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
                organism = (
                    str(organism).strip()
                    if pd.notna(organism) and organism
                    else "NA"
                )

                annotation = getattr(row, "Annotation", None)
                annotation = (
                    str(annotation).strip()
                    if pd.notna(annotation) and annotation
                    else "NA"
                )

                gram = getattr(row, "Gram_status", None)
                gram = (
                    str(gram).strip()
                    if pd.notna(gram) and gram
                    else "NA"
                )

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

        print(
            f"Saved {prefix} FASTA batch "
            f"{batch_idx + 1}/{n_batches}: {batch_file}"
        )


def export_pathogen_fastas(pathogen_subset: pd.DataFrame) -> None:
    gram_pos_df = pathogen_subset[pathogen_subset["Gram_status"] == "positive"]
    gram_neg_df = pathogen_subset[pathogen_subset["Gram_status"] == "negative"]

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


def export_iedb_source_fastas(match_df: pd.DataFrame) -> None:
    iedb_fasta_out = DATA_DIR / "intermediate/matched_iedb_source_proteins.fasta"
    unresolved_out = DATA_DIR / "intermediate/unresolved_iedb_protein_ids.csv"

    print("Preparing matched IEDB source proteins...")

    matched_iedb_ids = (
        match_df["IEDB_Protein_ID"]
        .map(clean_id)
        .dropna()
        .unique()
        .tolist()
    )

    print(f"Unique matched IEDB source protein IDs: {len(matched_iedb_ids)}")

    resolved = 0
    unresolved = []

    with open(iedb_fasta_out, "w", encoding="utf-8") as out:
        for protein_id in tqdm(
            matched_iedb_ids,
            desc="Fetching IEDB proteins",
            unit="protein"
        ):
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
        save_csv(pd.DataFrame(unresolved), unresolved_out)
        print(f"Saved unresolved IDs to: {unresolved_out}")


def main() -> None:
    match_df, protein_meta = load_input_data()

    pathogen_subset = get_matched_pathogen_subset(
        match_df=match_df,
        protein_meta=protein_meta
    )

    pathogen_subset = add_gram_status(pathogen_subset)

    export_pathogen_fastas(pathogen_subset)

    export_iedb_source_fastas(match_df)

    print("Done.")




if __name__ == "__main__":
    main()