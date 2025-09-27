"""
callbacks.py
--------------
All Dash callbacks for:
- Chemical drill-down (summary + tabs: Concentration, pH, Co-occurrence)
- Protein drill-down table + AA composition
"""

import dash
from dash import html, dcc
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
from collections import Counter

from app.data_utils import parse_ph_series, concentration_series_mm_and_pct
from app.figures import (
    make_hist_with_kde_binwidth,
    cooccurrence_heatmap_and_topbar,
    empty_fig,
    aa_composition_bar,   # ✅ dark version comes from figures.py
)


# -----------------------------
# Amino acid utilities
# -----------------------------
def compute_sequence_length(fasta: str) -> int:
    """Return sequence length from FASTA string."""
    if not fasta or not isinstance(fasta, str):
        return 0
    return len("".join(fasta.split()))


def compute_aa_composition(fasta: str) -> dict:
    """Return normalized amino acid composition (20 standard AAs)."""
    aas = "ACDEFGHIKLMNPQRSTVWY"
    seq = "".join(fasta.split()).upper()
    counts = Counter(seq)
    total = sum(counts.values())
    if total == 0:
        return {aa: 0 for aa in aas}
    return {aa: counts.get(aa, 0) / total for aa in aas}


# -----------------------------
# Helpers
# -----------------------------
def quant_summary_with_counts(series: pd.Series) -> tuple[str, int]:
    """Return formatted summary and count for 5–95% range + median."""
    s = pd.to_numeric(series, errors="coerce").dropna()
    if s.empty:
        return "N/A", 0
    q05 = s.quantile(0.05)
    q50 = s.quantile(0.50)
    q95 = s.quantile(0.95)
    return f"{q05:.2f} – {q95:.2f} (median {q50:.2f})", len(s)


def summary_block(d: pd.DataFrame, df_all: pd.DataFrame, selected: str) -> html.Div:
    """Return HTML summary block with counts, CID, sequence length, and concentration ranges."""
    total_count = len(d)
    total_all = max(len(df_all), 1)
    percent = 100 * total_count / total_all
    unique_proteins = d["Protein_ID"].nunique()
    cid_vals = d["CID"].dropna().unique()
    cid_str = str(cid_vals[0]) if len(cid_vals) else "N/A"
    pubchem_link = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid_str}" if cid_str != "N/A" else None

    # --- Average sequence length ---
    if "FASTA_Sequence" in d.columns:
        lengths = d["FASTA_Sequence"].dropna().apply(len)
        avg_len = lengths.mean() if not lengths.empty else None
    else:
        avg_len = None

    # Concentration summaries
    mm_series, pct_series = concentration_series_mm_and_pct(d)
    mm_txt, mm_n = quant_summary_with_counts(mm_series)
    pct_txt, pct_n = quant_summary_with_counts(pct_series)

    bold_mm = mm_n >= pct_n and mm_n > 0
    bold_pct = pct_n >= mm_n and pct_n > 0

    mm_line = html.Div([
        html.Span("Typical mM (5–95%): "),
        html.Span(mm_txt, style={"fontWeight": "bold"} if bold_mm else {}),
        html.Span(f"  — n={mm_n}", className="text-muted ms-1"),
    ])

    pct_line = html.Div([
        html.Span("Typical % (5–95%): "),
        html.Span(pct_txt, style={"fontWeight": "bold"} if bold_pct else {}),
        html.Span(f"  — n={pct_n}", className="text-muted ms-1"),
    ])

    return html.Div(
        [
            html.Div(f"Occurrences: {total_count}  ({percent:.2f}% of all conditions)"),
            html.Div(f"Unique proteins: {unique_proteins}"),
            html.Div(
                ["CID: ", html.A(cid_str, href=pubchem_link, target="_blank")]
                if pubchem_link else [f"CID: {cid_str}"]
            ),
            html.Div(
                f"Average sequence length: {avg_len:.0f} aa"
                if avg_len is not None else "Average sequence length: N/A"
            ),
            html.Hr(className="my-2"),
            mm_line,
            pct_line,
        ],
        className="small",
    )


# -----------------------------
# Tabs content
# -----------------------------
def build_concentration_tab(d, selected, binwidth_mm, binwidth_pct, focus_iqr, logy):
    binwidth_mm = max(0.01, min(float(binwidth_mm), 10.0))
    binwidth_pct = max(0.005, min(float(binwidth_pct), 3.0))

    mm_series, pct_series = concentration_series_mm_and_pct(d)

    fig_mm = make_hist_with_kde_binwidth(
        mm_series, f"{selected} – Concentration (mM)", "mM",
        bin_width=binwidth_mm, focus_iqr=focus_iqr, logy=logy,
    )
    fig_pct = make_hist_with_kde_binwidth(
        pct_series, f"{selected} – Concentration (%)", "%",
        bin_width=binwidth_pct, focus_iqr=focus_iqr, logy=logy,
    )

    return dbc.Row(
        [
            dbc.Col(dbc.Card(dbc.CardBody([dcc.Graph(figure=fig_mm)]), class_name="mb-3"), width=6),
            dbc.Col(dbc.Card(dbc.CardBody([dcc.Graph(figure=fig_pct)]), class_name="mb-3"), width=6),
        ]
    )


def build_ph_tab(d, selected, focus_iqr):
    ph_numeric = parse_ph_series(d["pH"])
    s = pd.to_numeric(ph_numeric, errors="coerce").dropna()

    if s.empty:
        fig_ph = empty_fig("pH distribution")
    else:
        bw = 0.25  # fixed bin size
        hist = go.Histogram(
            x=s,
            xbins=dict(size=bw),
            marker=dict(line=dict(width=0)),
            opacity=0.85,
            name="Counts",
            marker_color="lightskyblue",
        )
        fig_ph = go.Figure(hist)
        fig_ph.update_layout(
            template="plotly_dark",  # ✅ match theme
            title=f"pH distribution for {selected}",
            xaxis_title="pH",
            yaxis_title="Count",
            margin=dict(l=40, r=20, t=50, b=40),
        )
        if focus_iqr and len(s) >= 3:
            q1, q3 = s.quantile(0.25), s.quantile(0.75)
            lo, hi = q1 - 1.5 * (q3 - q1), q3 + 1.5 * (q3 - q1)
            if lo < hi:
                fig_ph.update_xaxes(range=[max(s.min(), lo), min(s.max(), hi)])
        fig_ph.update_xaxes(type="linear", dtick=bw)

    return dbc.Row([dbc.Col(dbc.Card(dbc.CardBody([dcc.Graph(figure=fig_ph)])), width=12)])


def build_co_tab(df_all, d, selected):
    heatmap_fig, bar_fig = cooccurrence_heatmap_and_topbar(df_all, d, selected, top_k=15)
    return dbc.Row([
        dbc.Col(dbc.Card(dbc.CardBody([dcc.Graph(figure=heatmap_fig)])), width=12),
        dbc.Col(dbc.Card(dbc.CardBody([dcc.Graph(figure=bar_fig)])), width=12),
    ])


# -----------------------------
# Register callbacks
# -----------------------------
def register_callbacks(app, df):
    """Register all Dash callbacks."""

    @app.callback(
        dash.Output("chem-summary", "children"),
        dash.Output("tab-content", "children"),
        dash.Output("binwidth-mm-label", "children"),
        dash.Output("binwidth-pct-label", "children"),
        dash.Input("chem-dropdown", "value"),
        dash.Input("tabs", "value"),
        dash.Input("show-all", "value"),
        dash.Input("binwidth-mm", "value"),
        dash.Input("binwidth-pct", "value"),
    )
    def update_drilldown(selected, tab, show_all,
                         binwidth_mm, binwidth_pct):
        focus_iqr = "all" not in (show_all or [])

        bw_mm_clamped = max(0.01, min(float(binwidth_mm or 1.0), 10.0))
        bw_pct_clamped = max(0.005, min(float(binwidth_pct or 0.25), 3.0))
        mm_label = f"Bin width: {bw_mm_clamped:.2f} mM (cap: 10 mM)"
        pct_label = f"Bin width: {bw_pct_clamped:.2f} % (cap: 3 %)"

        if not selected:
            return (
                html.Div("⬆ Select a chemical to see details.", className="text-muted"),
                html.Div("No chemical selected.", className="text-muted"),
                mm_label, pct_label,
            )

        d = df[df["Standardized_Precipitate"] == selected]
        summary = summary_block(d, df, selected)

        if tab == "tab-conc":
            content = build_concentration_tab(d, selected, bw_mm_clamped, bw_pct_clamped, focus_iqr, logy=False)
        elif tab == "tab-ph":
            content = build_ph_tab(d, selected, focus_iqr)
        elif tab == "tab-co":
            content = build_co_tab(df, d, selected)
        else:
            content = html.Div("Select a tab.", className="text-muted")

        return summary, content, mm_label, pct_label

    @app.callback(
        dash.Output("protein-table", "children"),
        dash.Input("protein-submit", "n_clicks"),
        dash.State("protein-id-input", "value"),
    )
    def show_protein_conditions(n_clicks, protein_id):
        if not n_clicks or not protein_id:
            return html.Div("⬆ Enter a Protein_ID above and click Show Conditions.", className="text-muted")

        protein_id = protein_id.strip()
        d = df[df["Protein_ID"].astype(str).str.upper() == protein_id.upper()]
        if d.empty:
            return html.Div(f"No conditions found for Protein_ID '{protein_id}'.", className="text-danger")

        # Table of conditions
        table_df = d[["Standardized_Precipitate", "Concentration", "pH"]].copy()
        table_df = table_df.rename(columns={"Standardized_Precipitate": "Chemical"}).sort_values(
            ["Chemical", "Concentration"]
        )

        # Sequence stats
        seq = "".join(d["FASTA_Sequence"].dropna().astype(str))
        length = compute_sequence_length(seq)
        composition = compute_aa_composition(seq)
        aa_fig = aa_composition_bar(composition, f"Amino Acid Composition (Length={length})")

        return html.Div([
            dbc.Table.from_dataframe(table_df, striped=True, bordered=True, hover=True, responsive=True),
            dcc.Graph(figure=aa_fig)
        ])
