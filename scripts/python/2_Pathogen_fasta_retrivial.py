#%%
import os
import requests

# Input and output paths
proteome_id_file = "../Data/proteome_ids.txt"
output_dir = "../Data/proteome_fastas"

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

base_url = "https://rest.uniprot.org/uniprotkb/stream"

with open(proteome_id_file) as f:
    for line in f:
        proteome_id = line.strip()
        if not proteome_id:
            continue

        params = {
            "format": "fasta",
            "query": f"(proteome:{proteome_id})"
        }

        output_path = os.path.join(output_dir, f"{proteome_id}.fasta")

        response = requests.get(base_url, params=params, timeout=60)

        if response.status_code == 200 and response.text.strip():
            with open(output_path, "w") as out:
                out.write(response.text)
            print(f"Downloaded {proteome_id}")
        else:
            print(f"Failed to download {proteome_id} (status {response.status_code})")

# %%
