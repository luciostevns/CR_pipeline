import time
from io import StringIO

import requests
from Bio import SeqIO
from tqdm import tqdm

from crossreactivity.io import DATA_DIR, load_csv


BASE_URL = "https://rest.uniprot.org/uniprotkb/stream"


def load_proteome_ids() -> list[str]:
    proteome_ids_path = DATA_DIR / "intermediate/proteome_ids.csv"

    return (
        load_csv(proteome_ids_path)
        .iloc[:, 0]
        .dropna()
        .tolist()
    )


def fetch_proteome_fasta(proteome_id: str) -> str | None:
    """
    Fetch the FASTA sequence for a given proteome ID from the UniProt REST API.
    """

    params = {
        "format": "fasta",
        "query": f"(proteome:{proteome_id})",
    }

    try:
        response = requests.get(BASE_URL, params=params, timeout=300)

        if response.status_code != 200 or not response.text.strip():
            print(f"\nFailed: {proteome_id}")
            return None

        return response.text

    except Exception as e:
        print(f"\nError with {proteome_id}: {e}")
        return None


def write_proteome_records(fasta_text: str, proteome_id: str, output_handle) -> None:
    fasta_io = StringIO(fasta_text)

    for record in SeqIO.parse(fasta_io, "fasta"):
        record.description = f"{record.description} PROTEOME={proteome_id}"
        SeqIO.write(record, output_handle, "fasta")


def main() -> None:
    proteome_ids = load_proteome_ids()
    output_file = DATA_DIR / "raw/all_proteomes.fasta"

    print(f"Downloading {len(proteome_ids)} proteomes...")

    with open(output_file, "w", encoding="utf-8") as out:
        for proteome_id in tqdm(
            proteome_ids,
            desc="Downloading proteomes",
            unit="proteome"
        ):
            fasta_text = fetch_proteome_fasta(proteome_id)

            if fasta_text is None:
                continue

            write_proteome_records(
                fasta_text=fasta_text,
                proteome_id=proteome_id,
                output_handle=out
            )

            time.sleep(0.5)

    print(f"Saved: {output_file}")
    print("Done")


if __name__ == "__main__":
    main()