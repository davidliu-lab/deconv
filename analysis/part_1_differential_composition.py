# %%
# imports, ipython magics, logging
import logging

import duckdb
import plotly.express as px
from IPython.core.interactiveshell import InteractiveShell

import helpers
from helpers.deg_analysis import (
    classifier_metrics,
    displaying_tables,
    plotting_curves,
    plotting_utils,
    plotting_volcanos,
)
from helpers.deg_analysis.loading_results import (
    get_arrow_dataset_for_deg_analysis_results,
)
from helpers.deg_analysis.postprocessing_gene_stats_fields import add_more_pval_fields

if InteractiveShell.initialized():
    ipython = InteractiveShell.instance()
    ipython.run_line_magic("load_ext", "autoreload")
    ipython.run_line_magic("autoreload", "2")
    ipython.run_line_magic("sql", "duckdb:///:default:")

logging.basicConfig(level="INFO", format=helpers.logging.FORMAT)
logging.getLogger("duckdb").setLevel("INFO")
logging.getLogger("helpers").setLevel("DEBUG")

logger = logging.getLogger(__name__)
logger.setLevel("DEBUG")

# %%
# make arrow dataset
deg_analysis_results = get_arrow_dataset_for_deg_analysis_results(
    # "gs://liulab/differential_composition_and_expression/copied/20230505_21h41m44s/deg_analysis/"
    "gs://liulab/differential_composition_and_expression/20230615_01h18m52s/deg_analysis/"
)

# %%
# query arrow dataset
query_text = """
SELECT
    origin,
    malignant_means,
    log2_fc,  -- needed by add_more_pval_fields
    run_id,
    gene_symbol,
    perturbed AND log2_fc != 0 AS perturbed,
    log2_fold_change,
    "pval",
    "-log10_pval",
    "pval_adjusted_bh",
    -1.0 * log10("pval_adjusted_bh") as "-log10_pval_adjusted_bh",
    "significant_bh_fdr=0.10",
FROM deg_analysis_results
WHERE
    log2_fc = 0.0
    AND origin = 'malignant_cibersortx'
    --AND malignant_means in ('None,None', '0.65,0.75', '0.6,0.8') --, '0.57,0.83', '0.55,0.85')
    --AND run_id in (0, 1)
;
"""
df_gene_stats = duckdb.sql(query_text).df()
df_gene_stats = add_more_pval_fields(df_gene_stats)

# %%
# compute fpr manually from gene_stats for alpha = 0.1
df_table = classifier_metrics.get_metrics_for_alphas(
    df_gene_stats,
    ["malignant_means", "log2_fc", "run_id"],
    [0.05, 0.1, 0.25],
)

# %%

displaying_tables.make_score_table_with_stddev(
    df_table["fp"],
    cmap="Blues",
    index="malignant_means",
    columns="alpha",
)


# %%
# compute curves, new style
groupby_cols = ["malignant_means", "run_id"]
# score_col = "-log10_pval_adjusted_bh_signed_directional"
score_col = "-log10_pval_signed_directional"
df_curves = classifier_metrics.get_curves_with_all_pvals(
    df_gene_stats, groupby_cols, "perturbed", score_col
)


# %%
# volcano plots
df_to_plot = (
    df_gene_stats
    # .query("malignant_means in ('0.55,0.85', '0.6,0.8', '0.65,0.75', 'None,None')")
    .groupby(["malignant_means", "gene_symbol"])
    .median()
    .reset_index()
)
print(df_to_plot.groupby("malignant_means").size())
fig = px.scatter(
    df_to_plot,
    x="log2_fold_change",
    y="-log10_pval_adjusted_bh",
    hover_name="gene_symbol",
    hover_data=["-log10_pval"],
    facet_col="malignant_means",
    facet_col_wrap=4,
    title="Volcano plots",
)
fig = plotting_volcanos.format_volcano_figure(fig, marker_size=2)
fig.update_layout(width=1000, height=800)
# fig.show(renderer="png", scale=2)
fig


# %%
# plot curves
fig = px.line(
    df_curves.reset_index(),
    x="-log10_pval_adjusted_bh_signed_directional",
    y="fpr",
    facet_col="malignant_means",
    color="run_id",
    hover_data=["fp", "fpr"],
)
plotting_utils.add_fdr_lines(fig, alpha=0.1)
plotting_utils.remove_excess_facet_axis_titles(fig)


# %%
# FPR curves
fig = px.line(
    df_curves_signed_directional.reset_index(),
    x="-log10_pval_adjusted_bh_signed_directional",
    y="fpr",
    color="run_id",
    markers=True,
    facet_col="malignant_means",
)
fig = plotting_curves.format_metric_by_threshold_figure(fig)
fig.show(renderer="png", scale=2)


# %%
# PPV curves for different levels of differential composition
fig = plotting_curves.plot_metric_by_threshold(
    df_curves,
    score_column="-log10_pval_adjusted_bh_signed_directional",
    metric_column="precision",
    facet_col="malignant_means",
    facet_row=None,
)
# 1000 x 500
fig.update_layout(width=1000, height=500)
# png
fig.show(renderer="png", scale=2)


# %%
# volcano plots for different levels of differential composition

# from helpers.deg_analysis import plotting_volcanos

# fig = plotting_volcanos.make_volcano_grid_scatter(
#     df_gene_stats,
#     aggregate_by=["origin", "malignant_means", "log2_fc", "gene_symbol", "perturbed"],
#     pval_col="-log10_pval_adjusted_bh",
#     perturbed_col="perturbed",
#     facet_col="malignant_means",
#     facet_row=None,
#     marker_color="perturbed",
# )
# # 1000 x 500
# fig.update_layout(width=1000, height=500)
# # png
# fig.show(renderer="png", scale=2)
