from pathlib import Path
import pandas as pd
import requests

# Project paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"


def load_excel(path, **kwargs):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"{path} does not exist")

    # Try IEDB-style double header first, using calamine
    try:
        df = pd.read_excel(path, header=[0, 1], engine="calamine", **kwargs)

        if isinstance(df.columns, pd.MultiIndex):
            df.columns = [
                " - ".join([str(c).strip() for c in col if pd.notna(c)])
                for col in df.columns
            ]
        else:
            df.columns = df.columns.str.strip()

        return df

    except Exception:
        # Fallback to normal single-header Excel
        df = pd.read_excel(path, header=0, engine="calamine", **kwargs)
        df.columns = df.columns.str.strip()
        return df
from pathlib import Path
import pandas as pd


def load_csv(path, **kwargs):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"{path} does not exist")

    return pd.read_csv(
        path,
        sep=",",
        dtype=str,
        **kwargs
    )


def load_tsv(path, **kwargs):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"{path} does not exist")

    return pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        **kwargs
    )

def load_tsv_robust(path, **kwargs):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"{path} does not exist")

    return pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        engine="python",       # slower but tolerant
        on_bad_lines="skip",   # skip broken rows
        quoting=3,             # ignore quotes entirely
        **kwargs
    )

def save_csv(df, path, **kwargs):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, **kwargs)

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