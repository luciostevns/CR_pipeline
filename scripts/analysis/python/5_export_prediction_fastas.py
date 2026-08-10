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
    protein_sequences_file = DATA_DIR / "intermediate/protein_sequences.csv"

    print("Loading data files...")
    match_df = load_csv(match_file)
    protein_sequences = load_csv(protein_sequences_file)

    return match_df, protein_sequences


def prepare_pathogen_exports(
    match_df: pd.DataFrame,
    protein_sequences: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    match_cols = {"Pathogen_Protein_ID", "Pathogen_Organism"}
    sequence_cols = {"Protein_ID", "Sequence", "Annotation"}

    for name, df, required in [
        ("iedb_match_regions_long.csv", match_df, match_cols),
        ("protein_sequences.csv", protein_sequences, sequence_cols),
    ]:
        missing = required - set(df.columns)
        if missing:
            raise ValueError(f"{name} is missing columns: {sorted(missing)}")

    match_df = match_df.copy()
    protein_sequences = protein_sequences.copy()
    match_df["Protein_ID"] = match_df["Pathogen_Protein_ID"].map(clean_id)
    protein_sequences["Protein_ID"] = protein_sequences["Protein_ID"].map(clean_id)

    if match_df["Protein_ID"].isna().any():
        raise ValueError("Matched pathogen rows contain missing Protein_IDs.")
    if protein_sequences["Protein_ID"].isna().any():
        raise ValueError("protein_sequences.csv contains missing Protein_IDs.")
    if protein_sequences["Protein_ID"].duplicated().any():
        examples = protein_sequences.loc[
            protein_sequences["Protein_ID"].duplicated(keep=False), "Protein_ID"
        ].unique()[:10]
        raise ValueError(
            "protein_sequences.csv must contain one row per Protein_ID. "
            f"Duplicate examples: {', '.join(examples)}"
        )

    matched_ids = match_df[["Protein_ID"]].drop_duplicates()
    netsurfp = matched_ids.merge(
        protein_sequences[["Protein_ID", "Sequence", "Annotation"]],
        on="Protein_ID",
        how="left",
        validate="one_to_one",
        indicator=True,
    )
    missing_sequences = netsurfp.loc[netsurfp["_merge"].eq("left_only"), "Protein_ID"]
    if not missing_sequences.empty:
        raise ValueError(
            "Matched proteins missing from protein_sequences.csv: "
            + ", ".join(missing_sequences.head(10))
        )
    netsurfp = netsurfp.drop(columns="_merge")
    blank_sequences = (
        netsurfp["Sequence"].isna()
        | netsurfp["Sequence"].astype(str).str.strip().eq("")
    )
    if blank_sequences.any():
        raise ValueError(
            "Matched proteins with missing sequences: "
            + ", ".join(netsurfp.loc[blank_sequences, "Protein_ID"].head(10))
        )

    gram_membership = match_df[["Protein_ID", "Pathogen_Organism"]].drop_duplicates()
    gram_membership["Pathogen_Organism"] = (
        gram_membership["Pathogen_Organism"]
        .astype("string")
        .str.strip()
        .replace("", pd.NA)
    )
    gram_membership["Gram_status"] = gram_membership["Pathogen_Organism"].map(
        GRAM_STATUS
    )
    unmapped = gram_membership[gram_membership["Gram_status"].isna()]
    if not unmapped.empty:
        species = unmapped["Pathogen_Organism"].dropna().drop_duplicates().sort_values()
        missing_organisms = unmapped["Pathogen_Organism"].isna().sum()
        details = species.to_list()
        if missing_organisms:
            details.append(f"<missing organism in {missing_organisms} row(s)>")
        raise ValueError("Missing GRAM_STATUS mapping for: " + ", ".join(details))

    gram_membership = gram_membership[["Protein_ID", "Gram_status"]].drop_duplicates()
    invalid_gram = set(gram_membership["Gram_status"]) - {"positive", "negative"}
    if invalid_gram:
        raise ValueError(f"Unexpected Gram status values: {sorted(invalid_gram)}")

    deeploc = gram_membership.merge(
        netsurfp,
        on="Protein_ID",
        how="left",
        validate="many_to_one",
    )
    suffix = deeploc["Gram_status"].map({"positive": "GP", "negative": "GN"})
    deeploc["Prediction_ID"] = deeploc["Protein_ID"] + "_" + suffix
    manifest = deeploc[["Prediction_ID", "Protein_ID", "Gram_status"]].copy()

    cross_gram = (
        gram_membership.groupby("Protein_ID")["Gram_status"]
        .nunique()
        .gt(1)
        .sum()
    )
    print(f"Unique matched pathogen proteins for NetSurfP: {len(netsurfp):,}")
    print("DeepLoc protein/Gram-status inputs:")
    print(gram_membership["Gram_status"].value_counts().to_string())
    print(f"Proteins represented in both Gram classes: {cross_gram:,}")

    return netsurfp, deeploc, manifest


def export_fasta_batches(
    df: pd.DataFrame,
    prefix: str,
    id_column: str,
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
                prediction_id = str(getattr(row, id_column)).strip()
                annotation = getattr(row, "Annotation", None)
                annotation = (
                    str(annotation).strip()
                    if pd.notna(annotation) and annotation
                    else "NA"
                )

                sequence = getattr(row, "Sequence", None)
                sequence = str(sequence).strip()
                header = f"{prediction_id} annotation={annotation}"

                write_fasta_record(out, header, sequence)

        print(
            f"Saved {prefix} FASTA batch "
            f"{batch_idx + 1}/{n_batches}: {batch_file}"
        )


def export_pathogen_fastas(
    netsurfp: pd.DataFrame,
    deeploc: pd.DataFrame,
) -> None:
    export_fasta_batches(
        netsurfp,
        prefix="matched_pathogen_proteins_netsurfp",
        id_column="Protein_ID",
    )

    export_fasta_batches(
        deeploc[deeploc["Gram_status"].eq("positive")],
        prefix="matched_pathogen_proteins_deeploc_gram_positive",
        id_column="Prediction_ID",
    )

    export_fasta_batches(
        deeploc[deeploc["Gram_status"].eq("negative")],
        prefix="matched_pathogen_proteins_deeploc_gram_negative",
        id_column="Prediction_ID",
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
    match_df, protein_sequences = load_input_data()
    netsurfp, deeploc, manifest = prepare_pathogen_exports(
        match_df=match_df,
        protein_sequences=protein_sequences,
    )
    export_pathogen_fastas(netsurfp, deeploc)
    manifest_out = DATA_DIR / "intermediate/pathogen_deeploc_prediction_manifest.csv"
    save_csv(manifest, manifest_out)
    print(f"Saved DeepLoc prediction manifest to: {manifest_out}")

    export_iedb_source_fastas(match_df)

    print("Done.")




if __name__ == "__main__":
    main()