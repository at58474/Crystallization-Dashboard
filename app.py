"""
app.py
-------
Main entry point. Runs Dash server, builds layout, and registers callbacks.
"""

import dash
import dash_bootstrap_components as dbc
from app.data_utils import load_data
from app.figures import make_top50_overview
from app.layout import make_layout
from app.callbacks import register_callbacks
import plotly.io as pio

pio.templates.default = "plotly_dark"

# Load data
df = load_data()

# Dash app
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.SLATE])
app.title = "Crystallization EDA"

# 👇 expose the Flask server for Render
server = app.server

# Top 50 overview + chem dropdown
fig_top = make_top50_overview(df)
chem_options = [{"label": c, "value": c}
                for c in sorted(df["Standardized_Precipitate"].dropna().unique())]

# Layout
app.layout = make_layout(fig_top, chem_options)

# Callbacks
register_callbacks(app, df)

# Local run
if __name__ == "__main__":
    app.run(host="0.0.0.0", port=8050, debug=True)
