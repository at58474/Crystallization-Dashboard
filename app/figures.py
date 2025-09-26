"""
figures.py
-----------
Figure helper functions for Crystallization EDA dashboard.
All figures use dark theme (plotly_dark) to match SLATE Bootstrap theme.
"""

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


# --------------------
# Empty placeholder fig
# --------------------
def empty_fig(title: str):
    fig = go.Figure()
    fig.update_layout(
        title=title,
        template="plotly_dark",
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        annotations=[dict(
            text="No data", x=0.5, y=0.5, showarrow=False,
            xref="paper", yref="paper", font=dict(color="lightgray")
        )],
        margin=dict(l=40, r=20, t=50, b=40),
    )
    return fig


# --------------------
# Lightweight KDE
# --------------------
def kde_curve(x: np.ndarray, points: int = 200):
    """Lightweight Gaussian KDE (no SciPy)."""
    x = np.asarray(x, dtype=float)
    x = x[~np.isnan(x)]
    if x.size < 2:
        return None, None
    std = np.std(x)
    if std == 0:
        return None, None
    n = len(x)
    h = 1.06 * std * (n ** (-1 / 5))
    grid = np.linspace(x.min(), x.max(), points)
    diff = grid[:, None] - x[None, :]
    kernel = np.exp(-0.5 * (diff / h) ** 2) / (np.sqrt(2 * np.pi) * h)
    density = kernel.mean(axis=1)
    return grid, density


# --------------------
# Histograms with KDE overlay
# --------------------
def make_hist_with_kde_binwidth(series, title, xaxis, bin_width, focus_iqr, logy):
    s = pd.to_numeric(series, errors="coerce").dropna()
    if s.empty:
        return empty_fig(title)

    bw = float(bin_width) if bin_width and bin_width > 0 else (1.0 if xaxis == "mM" else 0.25)

    hist = go.Histogram(
        x=s,
        xbins=dict(size=bw),
        marker=dict(line=dict(width=0)),
        opacity=0.85,
        name="Counts",
        marker_color="lightskyblue",
    )
    fig = go.Figure(hist)
    fig.update_layout(
        template="plotly_dark",
        title=title,
        xaxis_title=xaxis,
        yaxis_title="Count",
        margin=dict(l=40, r=20, t=50, b=40),
        bargap=0.02,
    )

    gx, gy = kde_curve(s.values)
    if gx is not None and gy is not None and np.all(np.isfinite(gy)):
        data_range = s.max() - s.min()
        approx_bins = max(int(round(data_range / bw)), 1)
        approx_bin_height = len(s) / approx_bins if approx_bins > 0 else len(s)
        fig.add_trace(go.Scatter(
            x=gx, y=gy * approx_bin_height, mode="lines",
            name="KDE", line=dict(width=2, color="cyan")
        ))

    if focus_iqr and len(s) >= 3:
        q1, q3 = s.quantile(0.25), s.quantile(0.75)
        iqr = q3 - q1
        lo, hi = q1 - 1.5 * iqr, q3 + 1.5 * iqr
        if lo < hi:
            fig.update_xaxes(range=[max(s.min(), lo), min(s.max(), hi)])

    if logy:
        fig.update_yaxes(type="log")

    return fig


# --------------------
# Co-occurrence figures
# --------------------
def cooccurrence_heatmap_and_topbar(df_all, d_sel, selected, top_k=15):
    co_df = df_all[
        (df_all["Protein_ID"].isin(d_sel["Protein_ID"]))
        & (df_all["Standardized_Precipitate"] != selected)
    ].copy()

    if co_df.empty:
        return empty_fig("Co-occurrence heatmap"), empty_fig("Top 10 co-occurring chemicals")

    counts = (
        co_df.groupby("Standardized_Precipitate")
        .size()
        .reset_index(name="count")
        .sort_values("count", ascending=False)
    )

    top10 = counts.head(10)
    bar_fig = px.bar(
        top10,
        x="Standardized_Precipitate",
        y="count",
        template="plotly_dark",
        title=f"Top 10 co-occurring chemicals with {selected}",
        color="count",
        color_continuous_scale="Viridis",
    )
    bar_fig.update_layout(xaxis_tickangle=45, margin=dict(l=40, r=20, t=50, b=120))

    top_chems = counts["Standardized_Precipitate"].head(min(top_k, len(counts))).tolist()
    mat_df = (
        co_df[co_df["Standardized_Precipitate"].isin(top_chems)]
        .drop_duplicates(subset=["Protein_ID", "Standardized_Precipitate"])
        .assign(val=1)
        .pivot_table(index="Protein_ID", columns="Standardized_Precipitate", values="val", fill_value=0)
        .reindex(columns=top_chems, fill_value=0)
    )

    if mat_df.empty:
        return empty_fig("Co-occurrence heatmap"), bar_fig

    X = mat_df.values
    co_mat = X.T @ X
    heatmap_fig = go.Figure(data=go.Heatmap(
        z=co_mat, x=top_chems, y=top_chems,
        colorscale="Viridis", colorbar=dict(title="Count")
    ))
    heatmap_fig.update_layout(
        title=f"Co-occurrence heatmap (subset with {selected})",
        template="plotly_dark",
        xaxis=dict(tickangle=45),
        margin=dict(l=60, r=20, t=50, b=120),
    )

    return heatmap_fig, bar_fig


# --------------------
# Top 50 overview
# --------------------
def make_top50_overview(df: pd.DataFrame):
    overall_n = len(df)
    chem_counts = (
        df.groupby("Standardized_Precipitate", dropna=True)
        .agg(count=("Standardized_Precipitate", "size"),
             unique_proteins=("Protein_ID", pd.Series.nunique))
        .reset_index()
        .sort_values("count", ascending=False)
    )
    chem_counts["percent"] = (chem_counts["count"] / max(overall_n, 1) * 100.0).round(2)
    top50 = chem_counts.head(50).copy()

    fig_top = px.bar(
        top50,
        x="Standardized_Precipitate",
        y="count",
        template="plotly_dark",
        title="Top 50 Standardized Chemicals by Frequency",
        color="count",
        color_continuous_scale="Viridis",
    )
    fig_top.update_traces(
        hovertemplate="<b>%{x}</b><br>"
                      "Count: %{y}<br>"
                      "Unique proteins: %{customdata[0]}<br>"
                      "% of all: %{customdata[1]}%",
        customdata=np.stack([top50["unique_proteins"], top50["percent"]], axis=-1),
    )
    fig_top.update_layout(
        xaxis_title="Chemical",
        yaxis_title="Count",
        xaxis_tickangle=45,
        height=520,
        margin=dict(l=40, r=40, t=60, b=140),
    )
    return fig_top


# --------------------
# Amino acid composition chart
# --------------------
def aa_composition_bar(composition: dict, title: str):
    df = pd.DataFrame({"AA": list(composition.keys()), "Fraction": list(composition.values())})
    fig = px.bar(
        df,
        x="AA",
        y="Fraction",
        title=title,
        template="plotly_dark",
        color="Fraction",
        color_continuous_scale="Viridis",
    )
    fig.update_layout(
        yaxis=dict(range=[0, 1]),
        margin=dict(l=40, r=20, t=50, b=40),
        xaxis_title="Amino Acid",
        yaxis_title="Fraction",
    )
    return fig
