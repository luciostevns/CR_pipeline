import requests

def extract_genus_species(scientific_name):
    """
    Extract genus/species as the first two words of the scientific name.
    """
    if not scientific_name or not str(scientific_name).strip():
        return None

    parts = str(scientific_name).strip().split()

    if len(parts) >= 2:
        return " ".join(parts[:2])
    elif len(parts) == 1:
        return parts[0]
    else:
        return None


def fetch_proteome_metadata(proteome_id, session=None, timeout=30):
    """
    Fetch proteome-level metadata from UniProt REST.

    Returns a dict with:
    - Proteome_ID
    - Scientific_name
    - Genus_species
    - Strain
    """
    if session is None:
        session = requests.Session()

    url = f"https://rest.uniprot.org/proteomes/{proteome_id}"

    result = {
        "Proteome_ID": proteome_id,
        "Scientific_name": None,
        "Genus_species": None,
        "Strain": None
    }

    try:
        response = session.get(url, timeout=timeout)
        response.raise_for_status()

        data = response.json()

        scientific_name = data.get("taxonomy", {}).get("scientificName", None)
        strain = data.get("strain", None)

        result["Scientific_name"] = scientific_name
        result["Genus_species"] = extract_genus_species(scientific_name)
        result["Strain"] = strain

    except Exception as e:
        print(f"Warning: failed to fetch metadata for {proteome_id}: {e}")

    return result

def fetch_uniprot_fasta(protein_id: str, timeout: int = 60) -> str | None:
    url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.fasta"
    try:
        r = requests.get(url, timeout=timeout)
        if r.status_code == 200 and r.text.startswith(">"):
            return r.text
    except Exception:
        pass
    return None