"""
main.py
--------
Main entry point. Runs Dash server, builds layout, and registers callbacks.
"""

import dash
import dash_bootstrap_components as dbc
import plotly.io as pio

from app.data_utils import run_query
from app.figures import make_top50_overview
from app.layout import make_layout
from app.callbacks import register_callbacks

pio.templates.default = "plotly_dark"

# -----------------------------
# Load initial data
# -----------------------------
TABLE = "conditions"

# Build Top 50 overview figure
sql_top50 = f"""
    SELECT Standardized_Precipitate, Protein_ID
    FROM {TABLE}
    WHERE Standardized_Precipitate IS NOT NULL
"""
df_top50 = run_query(sql_top50)
fig_top = make_top50_overview(df_top50)

# Build dropdown options directly from DB
sql_opts = f"""
    SELECT DISTINCT Standardized_Precipitate
    FROM {TABLE}
    WHERE Standardized_Precipitate IS NOT NULL
"""
chem_values = run_query(sql_opts)["Standardized_Precipitate"].sort_values().tolist()
chem_options = [{"label": c, "value": c} for c in chem_values]

# -----------------------------
# Dash app
# -----------------------------
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.SLATE])
app.title = "Crystallization EDA"

# Layout
app.layout = make_layout(fig_top, chem_options)

# Callbacks
register_callbacks(app)

# Expose server for Render
server = app.server

# -----------------------------
# Run local
# -----------------------------
if __name__ == "__main__":
    app.run(debug=True)
