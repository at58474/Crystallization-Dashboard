"""
main.py
-------
Entry point for the Crystallization EDA dashboard.
- Configures Dash
- Builds initial layout
- Registers callbacks
"""

import dash
import dash_bootstrap_components as dbc
import plotly.io as pio

from app.figures import make_top50_overview
from app.layout import make_layout
from app.callbacks import register_callbacks
from app.data_utils import get_db_path, fetch_top50_chemicals

pio.templates.default = "plotly_dark"

# Dash app
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.SLATE])
app.title = "Crystallization EDA"

# 👇 expose Flask server for Render
server = app.server

# Top 50 overview + chem dropdown (query directly from DB)
df_top50 = fetch_top50_chemicals()
fig_top = make_top50_overview(df_top50)
chem_options = [{"label": c, "value": c}
                for c in sorted(df_top50["Standardized_Precipitate"].dropna().unique())]

# Layout
app.layout = make_layout(fig_top, chem_options)

# Callbacks — they will query DB on demand
register_callbacks(app)

# Local run
if __name__ == "__main__":
    app.run(debug=True)
