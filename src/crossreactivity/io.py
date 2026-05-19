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

def clean_id(x: str) -> str | None:
    if pd.isna(x):
        return None
    x = str(x).strip()
    if not x or x.lower() == "nan":
        return None
    return x.upper()