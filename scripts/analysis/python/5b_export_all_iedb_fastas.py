"""Export full source-protein FASTAs for all records in wrangled_IEDB.csv.

Run in the same environment as 5_export_prediction_fastas.py. Existing
all-IEDB and matched-IEDB FASTAs are reused; missing proteins are fetched.
Outputs are separate from the matched-only prediction inputs.
"""

import time

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

from crossreactivity.io import DATA_DIR, load_csv, save_csv, clean_id
from crossreactivity.fasta import write_fasta_record, parse_fasta_text
from crossreactivity.uniprot import fetch_uniprot_fasta


BATCH_SIZE = 500
REQUEST_DELAY = 0.2


def accession(value):
    """Extract the accession from a plain ID, URL, or FASTA identifier."""
    if pd.isna(value) or not str(value).strip():
        return None
    token = str(value).strip().split()[0].rsplit("/", 1)[-1]
    token = token.split("|")[1] if "|" in token else token
    result = clean_id(token)
    return str(result).strip().upper() if pd.notna(result) and str(result).strip() else None


def read_cached_sequences(paths, requested_ids):
    sequences = {}
    for path in paths:
        if not path.exists():
            continue
        for record in SeqIO.parse(path, "fasta"):
            protein_id = accession(record.id)
            if protein_id not in requested_ids:
                continue
            sequence = str(record.seq)
            if not sequence:
                continue
            if protein_id in sequences and sequences[protein_id] != sequence:
                raise ValueError(f"Conflicting cached sequences for {protein_id}: {path}")
            sequences[protein_id] = sequence
    return sequences


def main():
    intermediate = DATA_DIR / "intermediate"
    intermediate.mkdir(parents=True, exist_ok=True)
    combined_path = intermediate / "all_iedb_source_proteins.fasta"
    batch_dir = intermediate / "all_iedb_prediction_fastas"
    batch_dir.mkdir(parents=True, exist_ok=True)

    iedb_data = load_csv(intermediate / "wrangled_IEDB.csv")
    if "Protein_ID" not in iedb_data:
        raise ValueError("wrangled_IEDB.csv is missing Protein_ID.")
    normalized_ids = iedb_data["Protein_ID"].map(accession)
    missing_ids = normalized_ids.isna()
    save_csv(
        iedb_data.loc[missing_ids].copy(),
        intermediate / "all_iedb_rows_missing_protein_id.csv",
    )
    protein_ids = sorted(normalized_ids.dropna().unique())
    if not protein_ids:
        raise ValueError("No source-protein IDs available for FASTA export.")

    sequences = read_cached_sequences(
        [combined_path, intermediate / "matched_iedb_source_proteins.fasta"],
        set(protein_ids),
    )
    print(f"IEDB records without a protein ID: {missing_ids.sum():,}")
    print(f"Unique source proteins requested: {len(protein_ids):,}")
    print(f"Source proteins available in existing FASTAs: {len(sequences):,}")

    # Rebuild the cache from reusable records, then save each successful
    # download immediately so an interrupted run can resume from this file.
    unresolved = []
    with combined_path.open("w", encoding="utf-8") as out:
        for protein_id in protein_ids:
            if protein_id in sequences:
                write_fasta_record(out, protein_id, sequences[protein_id])
        out.flush()

        pending = [protein_id for protein_id in protein_ids if protein_id not in sequences]
        for protein_id in tqdm(pending, desc="Fetching all IEDB source proteins", unit="protein"):
            try:
                fasta_text = fetch_uniprot_fasta(protein_id)
                parsed = parse_fasta_text(fasta_text) if fasta_text else None
                if parsed is None:
                    raise ValueError("No FASTA sequence returned")
                header, sequence = parsed
                if not sequence:
                    raise ValueError("Empty FASTA sequence returned")
                returned_id = accession(header.lstrip(">"))
                if returned_id != protein_id:
                    raise ValueError(
                        f"Returned accession {returned_id} differs from requested {protein_id}"
                    )
                sequences[protein_id] = sequence
                write_fasta_record(out, protein_id, sequence)
                out.flush()
            except Exception as exc:
                unresolved.append({"IEDB_Protein_ID": protein_id, "Error": str(exc)})
            finally:
                time.sleep(REQUEST_DELAY)

    # Always overwrite the unresolved report, including when every ID resolves.
    save_csv(
        pd.DataFrame(unresolved, columns=["IEDB_Protein_ID", "Error"]),
        intermediate / "all_iedb_unresolved_protein_ids.csv",
    )

    resolved_ids = [protein_id for protein_id in protein_ids if protein_id in sequences]
    manifest_rows = []
    for offset in range(0, len(resolved_ids), BATCH_SIZE):
        batch_ids = resolved_ids[offset:offset + BATCH_SIZE]
        batch_path = batch_dir / f"all_iedb_source_proteins_part_{offset // BATCH_SIZE + 1:03d}.fasta"
        with batch_path.open("w", encoding="utf-8") as out:
            for protein_id in batch_ids:
                write_fasta_record(out, protein_id, sequences[protein_id])
                manifest_rows.append({
                    "IEDB_Protein_ID": protein_id,
                    "Protein_Length": len(sequences[protein_id]),
                    "FASTA_Batch": str(batch_path),
                })

    manifest_path = intermediate / "all_iedb_prediction_manifest.csv"
    save_csv(pd.DataFrame(manifest_rows, columns=[
        "IEDB_Protein_ID", "Protein_Length", "FASTA_Batch"
    ]), manifest_path)
    print(f"Resolved source proteins: {len(resolved_ids):,}")
    print(f"Unresolved source proteins: {len(unresolved):,}")
    print(f"Combined FASTA: {combined_path}")
    print(f"Batch manifest: {manifest_path}")
    print("Use only batch files listed in the current manifest; older batches may remain.")


if __name__ == "__main__":
    main()
