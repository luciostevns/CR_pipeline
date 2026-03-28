#%%
import os
import requests
from pathlib import Path
from tqdm import tqdm
import re
from time import sleep


# Paths
proteome_id_file = "../Data/proteome_ids.txt"
fasta_folder = Path("../Data/proteome_fastas")
output_folder = Path("../Data/proteome_fastas_strain")
output_folder.mkdir(parents=True, exist_ok=True)

def fetch_strain_name(proteome_id):
    url = f"https://rest.uniprot.org/proteomes/{proteome_id}"
    response = requests.get(url)
    if response.status_code != 200:
        return None
    try:
        data = response.json()
        organism = data.get("taxonomy", {}).get("scientificName", "").strip()
        strain = data.get("strain", None)

        if not strain or strain.strip() == "":
            print(f"⚠️ No strain listed for {proteome_id}")
            return organism  # Return just the organism

        return f"{organism} {strain}"
    except Exception as e:
        print(f"Failed to parse JSON for {proteome_id}: {e}")
        return None


# Helper: update OS= field in FASTA header
def modify_fasta_headers(fasta_path, new_os_name, output_path):
    with open(fasta_path, "r") as infile, open(output_path, "w") as outfile:
        for line in infile:
            if line.startswith(">"):
                # Replace only the OS=... part up to (but not including) OX=
                line = re.sub(r"OS=.*?OX=", f"OS={new_os_name} OX=", line)
            outfile.write(line)

# Load proteome IDs
with open(proteome_id_file) as f:
    proteome_ids = [line.strip() for line in f if line.strip()]

# Main loop with tqdm
for proteome_id in tqdm(proteome_ids, desc="Processing proteomes"):
    fasta_path = fasta_folder / f"{proteome_id}.fasta"
    output_path = output_folder / f"{proteome_id}.fasta"

    if not fasta_path.exists():
        print(f"FASTA not found for {proteome_id}")
        continue

    strain_name = fetch_strain_name(proteome_id)
    if not strain_name:
        print(f"Strain not found for {proteome_id}")
        continue

    modify_fasta_headers(fasta_path, strain_name, output_path)
    sleep(0.3)  # gentle on UniProt's servers

print("✅ All headers updated.")


# %%
