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


def clean_id(x: str) -> str | None:
    if pd.isna(x):
        return None
    x = str(x).strip()
    if not x or x.lower() == "nan":
        return None
    return x.upper()

def fetch_uniprot_fasta(protein_id: str, timeout: int = 60) -> str | None:
    url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.fasta"
    try:
        r = requests.get(url, timeout=timeout)
        if r.status_code == 200 and r.text.startswith(">"):
            return r.text
    except Exception:
        pass
    return None

def parse_fasta_text(fasta_text: str) -> tuple[str, str] | None:
    if not fasta_text:
        return None

    lines = [line.strip() for line in fasta_text.splitlines() if line.strip()]
    if not lines or not lines[0].startswith(">"):
        return None

    header = lines[0][1:]
    sequence = "".join(lines[1:])

    if not sequence:
        return None

    return header, sequence


def fasta_wrap(seq: str, width: int = 60) -> str:
    seq = str(seq).strip()
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def write_fasta_record(handle, header: str, sequence: str) -> None:
    handle.write(f">{header}\n")
    handle.write(f"{fasta_wrap(sequence)}\n")




gram_status = {
    # Gram-negative
    "Acinetobacter baumannii": "negative",
    "Acinetobacter lactucae": "negative",
    "Acinetobacter nosocomialis": "negative",
    "Acinetobacter oleivorans": "negative",
    "Acinetobacter pittii": "negative",
    "Acinetobacter seifertii": "negative",
    "Aeromonas dhakensis": "negative",
    "Aeromonas hydrophila": "negative",
    "Aeromonas veronii": "negative",
    "Burkholderia aenigmatica": "negative",
    "Burkholderia ambifaria": "negative",
    "Burkholderia cenocepacia": "negative",
    "Burkholderia cepacia": "negative",
    "Burkholderia contaminans": "negative",
    "Burkholderia latens": "negative",
    "Burkholderia metallica": "negative",
    "Burkholderia multivorans": "negative",
    "Burkholderia orbicola": "negative",
    "Burkholderia pseudomultivorans": "negative",
    "Burkholderia vietnamiensis": "negative",
    "Campylobacter concisus": "negative",
    "Campylobacter gracilis": "negative",
    "Campylobacter jejuni": "negative",
    "Campylobacter lanienae": "negative",
    "Campylobacter lari": "negative",
    "Campylobacter peloridis": "negative",
    "Campylobacter rectus": "negative",
    "Campylobacter showae": "negative",
    "Campylobacter upsaliensis": "negative",
    "Campylobacter ureolyticus": "negative",
    "Campylobacter volucris": "negative",
    "Citrobacter amalonaticus": "negative",
    "Citrobacter farmeri": "negative",
    "Citrobacter freundii": "negative",
    "Citrobacter koseri": "negative",
    "Citrobacter murliniae": "negative",
    "Citrobacter sedlakii": "negative",
    "Citrobacter werkmanii": "negative",
    "Citrobacter youngae": "negative",
    "Cronobacter malonaticus": "negative",
    "Cronobacter sakazakii": "negative",
    "Edwardsiella tarda": "negative",
    "Elizabethkingia anophelis": "negative",
    "Elizabethkingia meningoseptica": "negative",
    "Enterobacter asburiae": "negative",
    "Enterobacter bugandensis": "negative",
    "Enterobacter chengduensis": "negative",
    "Enterobacter chuandaensis": "negative",
    "Enterobacter cloacae": "negative",
    "Enterobacter hormaechei": "negative",
    "Enterobacter intestinihominis": "negative",
    "Enterobacter kobei": "negative",
    "Enterobacter mori": "negative",
    "Enterobacter roggenkampii": "negative",
    "Enterobacter sichuanensis": "negative",
    "Escherichia albertii": "negative",
    "Escherichia coli": "negative",
    "Haemophilus influenzae": "negative",
    "Klebsiella electrica": "negative",
    "Klebsiella michiganensis": "negative",
    "Klebsiella oxytoca": "negative",
    "Klebsiella pneumoniae": "negative",
    "Klebsiella quasipneumoniae": "negative",
    "Klebsiella variicola": "negative",
    "Kluyvera ascorbata": "negative",
    "Kluyvera cryocrescens": "negative",
    "Kluyvera georgiana": "negative",
    "Kluyvera intermedia": "negative",
    "Legionella bozemanae": "negative",
    "Legionella pneumophila": "negative",
    "Morganella morganii": "negative",
    "Neisseria bacilliformis": "negative",
    "Neisseria cinerea": "negative",
    "Neisseria elongata": "negative",
    "Neisseria gonorrhoeae": "negative",
    "Neisseria lactamica": "negative",
    "Neisseria meningitidis": "negative",
    "Neisseria oralis": "negative",
    "Neisseria perflava": "negative",
    "Neisseria polysaccharea": "negative",
    "Neisseria subflava": "negative",
    "Neisseria weaveri": "negative",
    "Photobacterium damselae": "negative",
    "Pluralibacter gergoviae": "negative",
    "Providencia alcalifaciens": "negative",
    "Providencia rettgeri": "negative",
    "Providencia stuartii": "negative",
    "Pseudomonas aeruginosa": "negative",
    "Pseudomonas putida": "negative",
    "Raoultella planticola": "negative",
    "Salmonella diarizonae": "negative",
    "Salmonella dublin": "negative",
    "Salmonella enterica": "negative",
    "Salmonella enteritidis": "negative",
    "Salmonella hadar": "negative",
    "Salmonella newport": "negative",
    "Salmonella oranienberg": "negative",
    "Salmonella paratyphi": "negative",
    "Salmonella thompson": "negative",
    "Salmonella typhimurium": "negative",
    "Serratia marcescens": "negative",
    "Serratia odorifera": "negative",
    "Serratia proteamaculans": "negative",
    "Serratia rubidaea": "negative",
    "Shewanella algae": "negative",
    "Shigella boydii": "negative",
    "Shigella dysenteriae": "negative",
    "Shigella flexneri": "negative",
    "Shigella sonnei": "negative",
    "Stenotrophomonas maltophilia": "negative",
    "Vibrio cholerae": "negative",
    "Vibrio cidicii": "negative",
    "Vibrio fluvialis": "negative",
    "Vibrio mimicus": "negative",
    "Vibrio navarrensis": "negative",
    "Vibrio vulnificus": "negative",
    "Yersinia enterocolitica": "negative",

    # Gram-positive
    "Bacillus cereus": "positive",
    "Bacillus thuringiensis": "positive",
    "Clostridioides difficile": "positive",
    "Clostridium botulinum": "positive",
    "Clostridium innocuum": "positive",
    "Clostridium perfringens": "positive",
    "Corynebacterium aurimucosum": "positive",
    "Corynebacterium minutissimum": "positive",
    "Corynebacterium singulare": "positive",
    "Corynebacterium striatum": "positive",
    "Enterococcus faecalis": "positive",
    "Enterococcus faecium": "positive",
    "Intestinibacter bartlettii": "positive",
    "Listeria innocua": "positive",
    "Listeria monocytogenes": "positive",
    "Mycobacterium canetti": "positive",
    "Mycobacterium haemophilum": "positive",
    "Mycobacterium leprae": "positive",
    "Mycobacterium lepromatosis": "positive",
    "Mycobacterium tuberculosis": "positive",
    "Paraclostridium sordellii": "positive",
    "Peptacetobacter hiranonis": "positive",
    "Sarcina ventriculi": "positive",
    "Staphylococcus aureus": "positive",
    "Staphylococcus schweitzeri": "positive",
    "Streptococcus agalactiae": "positive",
    "Streptococcus equi": "positive",
    "Streptococcus halitosis": "positive",
    "Streptococcus infantis": "positive",
    "Streptococcus mitis": "positive",
    "Streptococcus mutans": "positive",
    "Streptococcus oralis": "positive",
    "Streptococcus peroris": "positive",
    "Streptococcus pneumoniae": "positive",
    "Streptococcus pseudopneumoniae": "positive",
    "Streptococcus pyogenes": "positive",
    "Terrisporobacter othiniensis": "positive",
}