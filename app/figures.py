"""
figures.py
------------
All figure-building utilities for the Crystallization EDA dashboard.
Now queries the DB directly (via run_query) instead of requiring
preloaded DataFrames.
"""

import plotly.express as px
import plotly.graph_objects as go
import pandas as pd

from app.data_utils import (
    run_query,
    parse_ph_series,
    concentration_series_mm_and_pct,
    compute_aa_composition,
    compute_sequence_length,
)

TABLE = "conditions"


# -----------------------------
# Utility empty figure
# -----------------------------
def empty_fig(title: str = "No data"):
    fig = go.Figure()
    fig.update_layout(
        title=title,
        template="plotly_dark",
        xaxis={"visible": False},
        yaxis={"visible": False},
        annotations=[{"text": title, "xref": "paper", "yref": "paper",
                      "showarrow": False, "font": {"color": "gray", "size": 16}}],
    )
    return fig


# -----------------------------
# Top 50 overview
# -----------------------------
def make_top50_overview():
    """
    Build top 50 precipitate overview bar chart.
    """
    sql = f"""
        SELECT Protein_ID, Standardized_Precipitate
        FROM {TABLE}
        WHERE Standardized_Precipitate IS NOT NULL
    """
    df = run_query(sql)

    if df.empty:
        return empty_fig("Top 50 Precipitates (no data)")

    agg = (
        df.groupby("Standardized_Precipitate")
          .agg(count=("Protein_ID", "nunique"))
          .reset_index()
          .sort_values("count", ascending=False)
          .head(50)
    )

    fig = px.bar(
        agg,
        x="Standardized_Precipitate",
        y="count",
        color="count",  # <--- add this back
        title="Top 50 Most Common Precipitates",
        template="plotly_dark",
        color_continuous_scale="Viridis",
    )
    fig.update_layout(
        xaxis_tickangle=45,
        margin=dict(l=40, r=20, t=50, b=120),
    )
    return fig


# -----------------------------
# Histograms with KDE
# -----------------------------
def make_hist_with_kde_binwidth(series, title, xlabel,
                                bin_width=1.0, focus_iqr=False, logy=False):
    s = pd.to_numeric(series, errors="coerce").dropna()
    if s.empty:
        return empty_fig(title)

    hist = go.Histogram(
        x=s,
        xbins=dict(size=bin_width),
        marker=dict(line=dict(width=0)),
        opacity=0.85,
        name="Counts",
    )

    fig = go.Figure(hist)
    fig.update_layout(
        template="plotly_dark",
        title=title,
        xaxis_title=xlabel,
        yaxis_title="Count",
        margin=dict(l=40, r=20, t=50, b=40),
    )

    if focus_iqr and len(s) >= 3:
        q1, q3 = s.quantile(0.25), s.quantile(0.75)
        lo, hi = q1 - 1.5*(q3-q1), q3 + 1.5*(q3-q1)
        fig.update_xaxes(range=[max(s.min(), lo), min(s.max(), hi)])

    if logy:
        fig.update_yaxes(type="log")

    return fig


# -----------------------------
# Co-occurrence heatmap + topbar
# -----------------------------
def cooccurrence_heatmap_and_topbar(df_all, df_selected, selected, top_k=15):
    """
    Given all Protein_ID–Precipitate pairs and selected subset,
    build a co-occurrence heatmap + top co-occurring bar chart.
    """
    if df_selected.empty or df_all.empty:
        return empty_fig("Co-occurrence heatmap"), empty_fig("Top co-occurring chemicals")

    proteins = df_selected["Protein_ID"].unique()
    co_df = df_all[df_all["Protein_ID"].isin(proteins)]

    counts = (
        co_df.groupby("Standardized_Precipitate")["Protein_ID"]
        .nunique()
        .reset_index(name="count")
        .sort_values("count", ascending=False)
    )

    top = counts.head(top_k)
    bar_fig = px.bar(
        top,
        x="Standardized_Precipitate",
        y="count",
        title=f"Top {top_k} Co-occurring Chemicals with {selected}",
        template="plotly_dark",
    )
    bar_fig.update_layout(xaxis_tickangle=45, margin=dict(l=40, r=20, t=50, b=120))

    # Heatmap: simple pivot table (chemicals × proteins)
    pivot = (
        co_df.pivot_table(
            index="Protein_ID",
            columns="Standardized_Precipitate",
            aggfunc="size",
            fill_value=0,
        )
    )
    heatmap_fig = px.imshow(
        pivot.T,
        aspect="auto",
        title=f"Co-occurrence Heatmap (proteins × chemicals) for {selected}",
        template="plotly_dark",
        color_continuous_scale="Viridis",
    )

    return heatmap_fig, bar_fig
