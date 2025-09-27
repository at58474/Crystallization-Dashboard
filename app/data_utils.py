"""
data_utils.py
----------------
Handles database loading and small preprocessing helpers for crystallization EDA.
"""

from pathlib import Path
import os
import sqlite3
import pandas as pd
import plotly.express as px

# Default: look for data/CrystallizationEDA.db relative to project root
ROOT_DIR = Path(__file__).resolve().parents[1]
DEFAULT_DB_PATH = ROOT_DIR / "data" / "CrystallizationEDA.db"
TABLE = "conditions"


def load_data(table: str = TABLE) -> pd.DataFrame:
    """
    Load conditions table from SQLite into a DataFrame.

    - Uses DB_PATH env var if set (Render deployment).
    - Falls back to DEFAULT_DB_PATH locally.
    - If no DB found, return empty DataFrame so app can still boot.
    """
    db_path = os.environ.get("DB_PATH", str(DEFAULT_DB_PATH))

    if not Path(db_path).exists():
        print(f"[WARN] Database not found at {db_path}. Starting with empty DataFrame.")
        # Minimal empty schema so the app doesn't break
        return pd.DataFrame(columns=["Protein_ID", "Standardized_Precipitate"])

    with sqlite3.connect(db_path) as conn:
        df = pd.read_sql(f"SELECT * FROM {table}", conn)
    return df


def parse_ph_series(series: pd.Series) -> pd.Series:
    """Extract numeric pH from strings like 'PH 7.5' or '7.5'."""
    def _one(s):
        if pd.isna(s):
            return None
        x = str(s).strip().upper()
        if x.startswith("PH"):
            x = x.replace("PH", "").strip()
        try:
            return float(x)
        except Exception:
            return None
    return series.apply(_one)


def concentration_series_mm_and_pct(d: pd.DataFrame):
    """
    Return (mm_series, pct_series):
     - mM = Concentration_Value where unit == 'mm'
       + Concentration_Converted for any rows whose unit is not 'mm' or '%'
     - %  = Concentration_Value where unit == '%'
    """
    unit = d["Concentration_Unit"].fillna("").astype(str).str.strip()
    unit_lower = unit.str.lower()

    mm_series = pd.concat(
        [
            d.loc[unit_lower == "mm", "Concentration_Value"],
            d.loc[~unit_lower.isin(["mm", "%"]), "Concentration_Converted"],
        ],
        ignore_index=True,
    )

    pct_series = d.loc[unit == "%", "Concentration_Value"]
    return mm_series, pct_series


def compute_sequence_length(fasta: str) -> int:
    if not fasta or not isinstance(fasta, str):
        return 0
    return len("".join(fasta.split()))


def compute_aa_composition(fasta: str) -> dict:
    from collections import Counter
    aas = "ACDEFGHIKLMNPQRSTVWY"  # standard 20 AA
    seq = "".join(fasta.split()).upper()
    counts = Counter(seq)
    total = sum(counts.values())
    if total == 0:
        return {aa: 0 for aa in aas}
    return {aa: counts.get(aa, 0) / total for aa in aas}


def aa_composition_bar(composition: dict, title: str):
    """
    Make a bar chart of amino acid composition with dark theme.
    """
    df = pd.DataFrame({
        "Amino Acid": list(composition.keys()),
        "Fraction": list(composition.values())
    })
    fig = px.bar(
        df,
        x="Amino Acid",
        y="Fraction",
        template="plotly_dark",   # 🔹 force dark background
        title=title,
        color="Fraction",
        color_continuous_scale="Viridis",
    )
    fig.update_layout(
        margin=dict(l=40, r=20, t=50, b=40),
        xaxis_title="Amino Acid",
        yaxis_title="Fraction",
    )
    return fig
