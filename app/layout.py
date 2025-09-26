"""
layout.py
-----------
Improved layout for Crystallization EDA dashboard:
- Theme tuned for SLATE (dark with contrast).
- Navbar header added.
- Protein section moved up (below Top 50).
- Chemical + View side by side above Drill-down.
- Drill-down spans full width.
- pH bin size fixed (0.25) → no slider.
- Removed log-scale Y axis toggle.
- Dark-themed dropdown styling (via CSS).
"""

import dash_bootstrap_components as dbc
from dash import dcc, html


def Card(children, title=None, dark=True):
    """Helper to wrap children in a styled Bootstrap Card for SLATE theme."""
    return dbc.Card(
        [
            dbc.CardHeader(title, className="fw-bold") if title else None,
            dbc.CardBody(children),
        ],
        class_name="mb-3 shadow-sm",
        color="dark" if dark else "light",
        inverse=dark,
        style={"borderRadius": "12px"},
    )


def make_layout(fig_top, chem_options):
    """Return the full Dash layout for the app."""
    return dbc.Container(
        [
            # ---------------- Navbar ----------------
            dbc.NavbarSimple(
                brand="Crystallization EDA Dashboard",
                color="primary",
                dark=True,
                fluid=True,
                className="mb-4 shadow",
            ),

            # ---------------- Overview ----------------
            dbc.Row(
                dbc.Col(Card(dcc.Graph(figure=fig_top, id="top-chemicals"), title="Top 50 Chemicals"), width=12),
                class_name="mb-2",
            ),

            # ---------------- Protein Drill-down (moved up) ----------------
            dbc.Row(
                dbc.Col(
                    Card(
                        [
                            html.H5("Protein Drill-Down", className="mb-3 text-info"),
                            html.Div(
                                [
                                    dbc.Input(
                                        id="protein-id-input",
                                        type="text",
                                        placeholder="Enter Protein_ID (e.g. 1AKH)",
                                        style={
                                            "maxWidth": "340px",
                                            "display": "inline-block",
                                            "marginRight": "10px",
                                            "backgroundColor": "#2c3e50",  # darker background
                                            "color": "white",              # white text
                                        },
                                    ),
                                    dbc.Button(
                                        "Show Conditions",
                                        id="protein-submit",
                                        n_clicks=0,
                                        color="info",
                                        className="shadow-sm",
                                    ),
                                ],
                                className="mb-3",
                            ),
                            html.Div(
                                id="protein-table",
                                style={"maxHeight": "520px", "overflowY": "auto"},
                            ),
                        ],
                        title="Protein",
                    ),
                    width=12,
                ),
                class_name="mt-4",
            ),

            # ---------------- Chemical + View side by side ----------------
            dbc.Row(
                [
                    dbc.Col(
                        Card(
                            [
                                html.Label("Select a chemical", className="fw-semibold"),
                                dcc.Dropdown(
                                    id="chem-dropdown",
                                    options=chem_options,
                                    placeholder="Choose a chemical...",
                                    persistence=True,
                                    style={
                                        "backgroundColor": "#2c3e50",
                                        "color": "white",
                                    },
                                    className="custom-dropdown",  # custom CSS target
                                ),
                                html.Hr(),
                                html.Div(id="chem-summary", className="small"),
                            ],
                            title="Chemical",
                        ),
                        width=6,
                    ),
                    dbc.Col(
                        Card(
                            [
                                html.Div("Concentration view", className="fw-semibold mb-2"),
                                dbc.Checklist(
                                    id="show-all",
                                    options=[{"label": "Show all data (disable IQR focus)", "value": "all"}],
                                    value=[],
                                    switch=True,
                                ),
                                html.Hr(),
                                html.Div("Bin width (mM)", className="small mb-1"),
                                dcc.Slider(
                                    id="binwidth-mm",
                                    min=0.05, max=10.0, step=0.5, value=1.0,
                                    tooltip={"placement": "bottom"},
                                ),
                                html.Div(
                                    id="binwidth-mm-label",
                                    className="text-muted small mt-1",
                                    children="Bin width: 1.00 mM (cap: 10 mM)",
                                ),
                                html.Div("Bin width (%)", className="small mt-3 mb-1"),
                                dcc.Slider(
                                    id="binwidth-pct",
                                    min=0.01, max=5.0, step=0.5, value=1.00,
                                    tooltip={"placement": "bottom"},
                                ),
                                html.Div(
                                    id="binwidth-pct-label",
                                    className="text-muted small mt-1",
                                    children="Bin width: 0.25 % (cap: 3 %)",
                                ),
                            ],
                            title="View",
                        ),
                        width=6,
                    ),
                ],
                class_name="mt-2",
            ),

            # ---------------- Drill-down full width ----------------
            dbc.Row(
                dbc.Col(
                    Card(
                        [
                            dcc.Tabs(
                                id="tabs",
                                value="tab-conc",
                                children=[
                                    dcc.Tab(label="Concentration", value="tab-conc", className="custom-tab"),
                                    dcc.Tab(label="pH", value="tab-ph", className="custom-tab"),
                                    dcc.Tab(label="Co-occurrence", value="tab-co", className="custom-tab"),
                                ],
                                className="custom-tabs",
                            ),
                            html.Div(id="tab-content", className="mt-3"),
                        ],
                        title="Drill-down",
                    ),
                    width=12,
                ),
                class_name="mt-3",
            ),
        ],
        fluid=True,
    )
