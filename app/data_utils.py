"""
data_utils.py
----------------
Handles database access and preprocessing helpers for crystallization EDA.
Now designed to query DB on demand instead of loading everything upfront.
"""

from pathlib import Path
import os
import sqlite3
import pandas as pd
import plotly.express as px


# -----------------------------
# DB Path Helpers
# -----------------------------
ROOT_DIR = Path(__file__).resolve().parents[1]
DEFAULT_DB_PATH = ROOT_DIR / "data" / "CrystallizationEDA.db"
TABLE = "conditions"


def get_db_path() -> str:
    """Return DB path (env var first, fallback to local)."""
    return os.environ.get("DB_PATH", str(DEFAULT_DB_PATH))


def run_query(sql: str, params: tuple = ()) -> pd.DataFrame:
    """Run a SQL query against the SQLite DB and return a DataFrame."""
    db_path = get_db_path()
    if not Path(db_path).exists():
        print(f"[WARN] Database not found at {db_path}. Returning empty DataFrame.")
        return pd.DataFrame()
    with sqlite3.connect(db_path) as conn:
        return pd.read_sql(sql, conn, params=params)


def fetch_top50_chemicals() -> pd.DataFrame:
    """Return counts of top 50 most common Standardized_Precipitates."""
    sql = f"""
        SELECT Standardized_Precipitate, COUNT(*) as count
        FROM {TABLE}
        WHERE Standardized_Precipitate IS NOT NULL
        GROUP BY Standardized_Precipitate
        ORDER BY count DESC
        LIMIT 50
    """
    return run_query(sql)


# -----------------------------
# Data Parsing Helpers
# -----------------------------
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
        template="plotly_dark",
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
