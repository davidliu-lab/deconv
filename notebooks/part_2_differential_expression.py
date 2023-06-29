# %%
# imports, logging, configuration
import logging
from turtle import width

import duckdb
from matplotlib import scale
from matplotlib.pyplot import sca
import numpy as np
import pandas as pd
import plotly.express as px
from IPython.core.interactiveshell import InteractiveShell
from google.cloud import bigquery

import helpers
from helpers.logging.configuring_logging import FORMAT, EasternTimeFormatter

logger = logging.getLogger(__name__)
logger.setLevel("DEBUG")
logging.getLogger("helpers").setLevel("INFO")
logging.getLogger("google.cloud").setLevel("INFO")
handler = logging.StreamHandler()
handler.setFormatter(EasternTimeFormatter(FORMAT, datefmt="%I:%M:%S %p"))
logging.getLogger().handlers = [handler]

if InteractiveShell.initialized():
    ipython = InteractiveShell.instance()
    ipython.run_line_magic("load_ext", "autoreload")
    ipython.run_line_magic("autoreload", "2")
    ipython.run_line_magic("sql", "duckdb:///:default:")

# %%
where_expression = """
    log2_fc != 0.0
    and malignant_means = 'None,None'
"""

# %%
# load gene stats

sql = f"""
SELECT *
FROM `deconv.deg_analysis_results`
WHERE {where_expression}
ORDER BY origin, malignant_means, log2_fc, run_id
;
"""
df_gene_stats = bigquery.Client().query(sql).to_dataframe(progress_bar_type="tqdm")


# %%
# metrics for different log2_fc values, alpha values
df_metrics = helpers.deg_analysis.classifier_metrics.get_metrics_for_alphas(
    df_gene_stats,
    groupbys=["origin", "log2_fc", "run_id"],
    score_col="neg_log10_pval_adjusted_bh_signed_directional",
)

# %% 
# analyzing bulk GEPs
helpers.deg_analysis.displaying_tables.make_score_table_with_stddev(
    df_metrics.query("origin == 'bulk'"),
    cmap="Blues",
    index="log2_fc",
    columns="alpha",
    vmin=0,
    vmax=5,
)

# %%
# analyzing cibersortx-inferred malignant-specific GEPs
helpers.deg_analysis.displaying_tables.make_score_table_with_stddev(
    df_metrics.query("origin == 'malignant_cibersortx'")["recall"],
    cmap="Blues",
    index="log2_fc",
    columns="alpha",
    vmin=0.0,
    vmax=5.0,
)

# %%
# ROC AUC box plot
fig = px.box(
    df_metrics.query("alpha == 0.05").reset_index(),
    x="log2_fc",
    y="roc_auc",
    color="origin",
    title="ROC AUC for different log2_fc values",
)
fig.update_yaxes(range=[0, 1])
fig.update_layout(legend=dict(yanchor="bottom", y=0.1, xanchor="center", x=0.5))
fig.show(renderer="png", scale=2)

# %%
# compute curves
groupby_cols = ["origin", "malignant_means", "log2_fc", "run_id"]
score_col = "neg_log10_pval_signed_directional"
df_curves = helpers.deg_analysis.classifier_metrics.get_curves_with_all_pvals(
    df_gene_stats, groupby_cols, "perturbed", score_col
)
df_curves

# %%
# plot ROC curves
fig = px.scatter(
    df_curves.query((
        "log2_fc in (-0.5, 0.5, 1.0)"
        # " and run_id < 3"
    )).reset_index(),
    x="fpr",
    y="tpr",
    facet_col="log2_fc",
    color="origin",
    trendline="lowess",
    trendline_options=dict(frac=0.1),
)
# set marker size to 1
fig.for_each_trace(lambda trace: trace.update(marker=dict(size=1)))
fig = helpers.deg_analysis.plotting_utils.remove_excess_facet_axis_titles(fig)
fig.update_xaxes(
    range=[0, 1],
    constrain="domain",
    dtick=0.2,
)
fig.update_yaxes(
    range=[0, 1],
    constrain="domain",
    dtick=0.2,
)
fig.update_layout(width=800, height=400)
fig.update_layout(legend=dict(yanchor="bottom", y=0.1, xanchor="right", x=0.99))
fig.update_layout(title="ROC curves for different log2_fc values")
fig.show(renderer="png", scale=2)
