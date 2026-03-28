from pathlib import Path
import pandas as pd


# =========================
# File Import and export helper functions
# =========================
def load_csv(path, **kwargs):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"{path} does not exist")
    return pd.read_csv(path, **kwargs)


def load_tsv(path, **kwargs):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"{path} does not exist")
    return pd.read_csv(path, sep="\t", **kwargs)


def save_csv(df, path, **kwargs):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, **kwargs)


# Validate that required columns are present
def require_columns(df: pd.DataFrame, required_cols: list[str], df_name: str = "DataFrame"):

    missing = [col for col in required_cols if col not in df.columns]

    if missing:
        raise ValueError(
            f"{df_name} is missing required columns: {missing}\n"
            f"Available columns: {list(df.columns)}"
        )