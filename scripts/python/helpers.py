from pathlib import Path
import pandas as pd



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

# Import and export helper functions
def load_csv(path, sep=",", **kwargs):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"{path} does not exist")
    return pd.read_csv(path, sep=sep, **kwargs)


def load_tsv(path, **kwargs):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"{path} does not exist")
    return pd.read_csv(path, sep="\t", **kwargs)
 

def save_csv(df, path, **kwargs):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, **kwargs)


# Validate if the correct cols are included
def require_columns(df, required_cols, df_name="DataFrame"):
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        raise ValueError(
            f"{df_name} is missing required columns: {missing}\n"
            f"Available columns: {list(df.columns)}"
        )