# %%
# imports, logging, configuration
import logging
from turtle import width

import duckdb
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
    log2_fc = 0.0
    AND malignant_means in (
        '0.55,0.85',
        '0.57,0.83',
        '0.6,0.8',
        '0.63,0.77',
        '0.65,0.75',
        '0.7,0.72',
        '0.71,0.71',
        '0.72,0.7',
        '0.75,0.65',
        '0.77,0.63',
        '0.8,0.6',
        '0.83,0.57',
        '0.85,0.55'
    )
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
# compute metrics for all runs, three alpha values
df_metrics = helpers.deg_analysis.classifier_metrics.get_metrics_for_alphas(
    df_gene_stats,
    groupbys=["origin", "malignant_means", "run_id"],
    score_col="neg_log10_pval_adjusted_bh",
)
# %%
# table of FP for all runs, three alpha values, bulk only
helpers.deg_analysis.displaying_tables.make_score_table_with_stddev(
    df_metrics.query("origin == 'bulk'")["fp"],
    cmap="Reds",
    index="malignant_means",
    columns="alpha",
    vmin=0,
    vmax=5000,
)
# %%
# table of FP for all runs, three alpha values, cibersortx only
helpers.deg_analysis.displaying_tables.make_score_table_with_stddev(
    df_metrics.query("origin == 'malignant_cibersortx'")["fp"],
    cmap="Reds",
    index="malignant_means",
    columns="alpha",
    vmin=0,
    vmax=5000,
)

# %%
# box plot of FP counts for all runs, three alpha values, bulk only
fig = px.box(
    df_metrics.query("alpha == 0.1").reset_index(),
    x="malignant_means",
    y="fp",
    color="origin",
)
fig.update_layout(
    width=800,
    height=500,
    title="False positives using FDR alpha=0.1 for different proportion scenarios",
    legend=dict(
        yanchor="top",
        y=0.95,
        xanchor="center",
        x=0.5,
    )
)
fig.show(renderer="png", scale=2)

# %%
# volcano plot of three malignant_means scenarios, cibersortx only
fig = px.scatter(
    df_gene_stats.query((
        # "origin == 'malignant_cibersortx'"
        # " and "
        "run_id == 1"
        " and malignant_means in ('0.57,0.83', '0.6,0.8', '0.63,0.77')"
    )),
    x="log2_fold_change",
    y="neg_log10_pval_adjusted_bh",
    facet_col="malignant_means",
    facet_row="origin"
)
fig = helpers.deg_analysis.plotting_volcanos.format_volcano_figure(fig, simplify_facet_titles=True, marker_size=3)
fig.update_xaxes(range=[-5, 5])
fig.update_yaxes(title_text="")
fig.update_layout(width=800, height=400)
fig.show(renderer="png", scale=2)


# %%
# volcano plot of origin == "malignant_cibersortx" and run_id == 0
fig = px.scatter(
    df_gene_stats.query("origin == 'malignant_cibersortx'"),
    x="log2_fold_change",
    y="neg_log10_pval_adjusted_bh",
    facet_col="malignant_means",
    # facet_col="run_id",
    facet_col_wrap=4,
    hover_name="gene_symbol",
)
fig = helpers.deg_analysis.plotting_volcanos.format_volcano_figure(fig, simplify_facet_titles=False)
fig.update_xaxes(range=[-4, 4])
fig.update_layout(
    width=800,
    height=1000,
    title="Data: CIBERSORTx-inferred malignant GEPs"
)
# remove y-axis label
# fig.update_yaxes(title_text="")
fig.show(renderer="png", scale=2)

# %%
# FP curve for all runs, by neg_log10_pval_adjusted_bh_signed_directional
df_curves = helpers.deg_analysis.classifier_metrics.get_curves_with_all_pvals(
     df_gene_stats.query((
        # "origin == 'malignant_cibersortx'"
        # " and "
        # "run_id == 1"
        # " and "
        "malignant_means in ('0.57,0.83', '0.6,0.8', '0.63,0.77')"
    )),
    groupby_cols=["origin", "malignant_means", "run_id"],
    score_col="neg_log10_pval_adjusted_bh",
)
fig = px.line(
    df_curves.reset_index(),
    x="neg_log10_pval_adjusted_bh",
    y="fp",
    color="run_id",
    # markers=True,
    facet_row="origin",
    facet_col="malignant_means",
)
fig = helpers.deg_analysis.plotting_curves.format_metric_by_threshold_figure(
    fig,
    add_fdr_lines=True,
    simplify_facet_titles=False,
)
fig.update_yaxes(range=[0, 5000])
fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
fig = helpers.deg_analysis.plotting_utils.remove_excess_facet_axis_titles(fig)
fig.show(renderer="png", scale=2)

# %%
