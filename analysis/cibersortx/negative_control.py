# %%
# imports, logging, configuration
import logging

import duckdb
import plotly.express as px

import helpers
from helpers.deg_analysis import plotting_curves, plotting_utils
from helpers.deg_analysis.computing_classifier_roc_precion_recall import (
    calculate_all_curves,
)
from helpers.deg_analysis.plotting_curves import plot_curves, plot_metric_by_threshold
from helpers.deg_analysis.plotting_utils import add_fdr_lines
from helpers.deg_analysis.plotting_volcanos_v1 import make_volcano_grid_scatter
from helpers.deg_analysis.postprocessing_gene_stats_fields import add_more_pval_fields
from helpers.running_cibersortx.loading_results import (
    get_arrow_dataset_for_deg_analysis_results,
)

logger = logging.getLogger(__name__)
logging.basicConfig(format=helpers.logging.LOGGING_FORMAT, level="INFO")

# %%
# load data
deg_analysis_results = get_arrow_dataset_for_deg_analysis_results(
    "gs://liulab/differential_composition_and_expression/copied/20230505_21h41m44s/deg_analysis/"
)

query_text = """
SELECT
    origin,
    malignant_means,
    log2_fc,
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
    malignant_means = 'None,None'
    --AND origin = 'malignant_cibersortx'
    AND log2_fc in (0.0)
;
"""
df_gene_stats = duckdb.sql(query_text).df()  # query with "duckdb:///:default:"
df_gene_stats = add_more_pval_fields(df_gene_stats)
df_curves = calculate_all_curves(
    df_gene_stats,
    score_col="-log10_pval_adjusted_bh_signed_directional",
    perturbed_col="perturbed",
)

# %%
# histogram of log2_fold_change for each origin
fig = px.histogram(
    df_gene_stats,
    x="log2_fold_change",
    facet_col="origin",
    nbins=100,
    histnorm="percent",
)
fig.update_xaxes(range=[-10, 10])
fig.update_yaxes(range=[0, 100])
fig.update_layout(width=800, height=400)
fig.show(renderer="png", scale=2)

# %%
# volcano plot: DGE for both origins
fig = make_volcano_grid_scatter(
    df=df_gene_stats,
    groupby_cols=["origin", "malignant_means", "log2_fc", "gene_symbol", "perturbed", "run_id"],
    perturbed_col="perturbed",
    facet_col="origin",
    pval_col="-log10_pval_adjusted_bh",
    simplify_facet_titles=False,
)
fig.update_xaxes(range=[-10, 10])
fig.update_yaxes(range=[0, 3])
fig.update_layout(
    title="Volcano plots DGE in simulated bulk RNA-seq and CIBERSORTx-inferred malignant GEPs"
)
fig.show(width=800, height=500, renderer="png", scale=2)

# %%
# volcano plot: DGE in simulated bulk data
fig = make_volcano_grid_scatter(
    # df=df_gene_stats.query("origin == 'bulk'"),
    df=df_gene_stats,
    groupby_cols=["origin", "malignant_means", "log2_fc", "gene_symbol", "perturbed", "run_id"],
    perturbed_col="perturbed",
    pval_col="-log10_pval_adjusted_bh",
    simplify_facet_titles=False,
)
fig.update_yaxes(range=[0, 3])
fig.update_layout(
    title="Volcano plot: DGE in simulated bulk data",
)
# fig.show(renderer="png", scale=2)

# %%
# volcano plot: DGE in CIBERSORTx-inferred malignant GEPs

fig = make_volcano_grid_scatter(
    df=df_gene_stats.query("origin == 'malignant_cibersortx'"),
    groupby_cols=["origin", "malignant_means", "log2_fc", "gene_symbol", "perturbed", "run_id"],
    perturbed_col="perturbed",
    pval_col="-log10_pval_adjusted_bh",
    simplify_facet_titles=False,
)
fig.update_yaxes(range=[0, 3])
fig.update_layout(
    title="Volcano plot: DGE in CIBERSORTx-inferred malignant GEPs",
)
# fig.show(renderer="png", scale=2)

# %%
# volcano plots of each run (DGE in CIBERSORTx-inferred malignant GEPs)

fig = make_volcano_grid_scatter(
    df=df_gene_stats.query("origin == 'malignant_cibersortx'"),
    groupby_cols=["origin", "malignant_means", "log2_fc", "gene_symbol", "perturbed", "run_id"],
    perturbed_col="perturbed",
    facet_col="run_id",
    pval_col="-log10_pval_adjusted_bh",
    simplify_facet_titles=False,
)
# make markers larger
fig.update_traces(marker=dict(size=5))
# remove axis limits for y-axis
# fig.update_yaxes(range=[None, None])
fig.update_layout(
    title="Volcano plot: DGE in CIBERSORTx-inferred malignant GEPs over several runs",
)
# remove all subplot x-axis titles
for i in range(1, 5):
    fig.layout[f"xaxis{i}"].title = None
# set x-axis title
fig.update_xaxes(title_text="log2_fold_change")
fig.update_layout(width=800, height=400)
fig.show(renderer="png", scale=2)

# %%
# skipping ROC curves.. not possible for no perturbation
# for each run (DGE in bulk and CIBERSORTx-inferred malignant GEPs)
"""
fig = px.line(
    df_curves.reset_index(),
    "fpr",
    "tpr",
    color="run_id",
    facet_col="origin",
    title="ROC curves for each run",
)
fig = plotting_curves.add_roc_dashed_lines(fig)
fig = plotting_curves.format_curves_fig(fig)
fig = plotting_utils.remove_variable_names_from_facet_axis_titles(fig)
fig.update_layout(width=800, height=500)
# fig.show(renderer="png", scale=2)
"""

# %%
# plot_curves(
#     df_curves.reset_index(),
#     facet_col="origin",
#     facet_row=None,
#     simplify_facet_titles=False,
# )

# %%
# plot FPR by BH-adjusted p-value
fig = plot_metric_by_threshold(
    df_curves=df_curves.query("origin == 'malignant_cibersortx'"),
    metric_column="fpr",
    score_column="-log10_pval_adjusted_bh_signed_directional",
    simplify_facet_titles=False,
)
fig.update_layout(
    title="FPR by BH-adjusted p-value",
)
fig.show(width=800, height=500, renderer="png", scale=2)

# %%
# box plot of p-values
fig = px.box(
    df_gene_stats.query("origin == 'malignant_cibersortx'"),
    x="origin",
    y="-log10_pval_adjusted_bh",
    color="perturbed",
    facet_col="log2_fc",
    facet_row="malignant_means",
)
# title: Distribution of BH-adjusted p-values for perturbed and non-perturbed genes
fig.update_layout(
    title="Distribution of BH-adjusted p-values for perturbed and non-perturbed genes",
)
fig = add_fdr_lines(
    fig=fig,
    alpha=0.1,
    log_scale=True,
    axis="y",
)
# legend overlay on top right
fig.update_layout(legend=dict(xanchor="right", yanchor="top", y=0.95, x=0.95))
fig.show(width=800, height=500, renderer="png", scale=5)

# %%
# volcano plots for all cell types (CIBERSORTx-inferred)

# create pyarrow dataset for GEP data

# compute DEG stats for all cell types

# plot volcano plots
