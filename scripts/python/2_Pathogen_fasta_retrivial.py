import requests
import time
from Bio import SeqIO
from io import StringIO
from tqdm import tqdm
from helpers import load_csv, DATA_DIR

proteome_ids = load_csv(DATA_DIR / "intermediate/proteome_ids.csv").iloc[:, 0].dropna().tolist()

output_file = DATA_DIR / "raw/all_proteomes.fasta"

base_url = "https://rest.uniprot.org/uniprotkb/stream"

print(f"Downloading {len(proteome_ids)} proteomes...")

with open(output_file, "w", encoding="utf-8") as out:

    for proteome_id in tqdm(proteome_ids, desc="Downloading proteomes", unit="proteome"):

        params = {
            "format": "fasta",
            "query": f"(proteome:{proteome_id})"
        }

        try:
            response = requests.get(base_url, params=params, timeout=300)

            if response.status_code != 200 or not response.text.strip():
                print(f"\nFailed: {proteome_id}")
                continue

            fasta_io = StringIO(response.text)

            for record in SeqIO.parse(fasta_io, "fasta"):
                record.description = f"{record.description} PROTEOME={proteome_id}"
                SeqIO.write(record, out, "fasta")

            time.sleep(0.5)

        except Exception as e:
            print(f"\nError with {proteome_id}: {e}")
            continue

print("Done")