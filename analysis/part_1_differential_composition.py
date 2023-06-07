# %%
# imports
import logging

import duckdb
import pandas as pd
import plotly.express as px

import helpers
from helpers.deg_analysis import displaying_tables, plotting_curves, plotting_volcanos
from helpers.deg_analysis.classifier_metrics import (
    calculate_all_curves,
    compute_scores,
)
from helpers.deg_analysis.postprocessing_gene_stats_fields import add_more_pval_fields
from helpers.running_cibersortx.loading_results import get_arrow_dataset_for_deg_analysis_results


# %%
# logging
logging.basicConfig(level="INFO", format=helpers.logging.FORMAT)
logging.getLogger("duckdb").setLevel("INFO")
logging.getLogger("helpers").setLevel("DEBUG")

logger = logging.getLogger(__name__)
logger.setLevel("DEBUG")

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
    --malignant_means in ('None,None', '0.65,0.75', '0.6,0.8')
    origin = 'malignant_cibersortx'
    AND log2_fc in (0.0,)
    AND run_id in (0, 1)
;
"""
df_gene_stats = duckdb.sql(query_text).df()  # query with "duckdb:///:default:"
df_gene_stats = add_more_pval_fields(df_gene_stats)
df_curves = calculate_all_curves(
    df_gene_stats,
    score_col="-log10_pval_adjusted_bh_signed_directional",
    perturbed_col="perturbed",
)
df_scores = compute_scores(df_gene_stats, perturbed_col="perturbed")


# %%
# volcano plots
import importlib

importlib.reload(plotting_volcanos)

plotting_volcanos.make_volcano_grid_scatter(
    df_gene_stats,
    groupby_cols=[
        "origin",
        "malignant_means",
        "log2_fc",
        "gene_symbol",
        "perturbed",
    ],
    pval_col="-log10_pval_adjusted_bh",
    hover_data=["-log10_pval"],
    perturbed_col="perturbed",
    facet_col="malignant_means",
    facet_row=None,
)


# %%
df_scores


# %%
# FP count at FDR alpha=0.1
def get_metrics_at_threshold(
    df: pd.DataFrame,
    groupby: list[str],
    score: str,
    threshold: float,
    metrics: list[str],
) -> pd.DataFrame:
    df = df.reset_index()
    df = df[df[score] >= threshold]
    dfg = df.groupby(groupby)
    result = dfg.apply(lambda x: x.sort_values(score, ascending=False)[metrics].iloc[-1])
    return result


df_metrics = get_metrics_at_threshold(
    df_curves,
    groupby=df_curves.index.names[:-1],
    score="-log10_pval_adjusted_bh_signed_directional",
    thresholdAq1=1.0,
    metrics=["fp", "tn", "fpr"],
)

df_metrics


# displaying_metric_tables.make_score_table_with_stddev(
#     fp_counts,
# )

# %%
# histograms of PPV (precision) for each experiment
px.box(
    df_scores,
    y="precision",
    x="malignant_means",
    # facet_col="malignant_means",
)


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
from helpers.deg_analysis import plotting_volcanos

fig = plotting_volcanos.make_volcano_grid_scatter(
    df_gene_stats,
    groupby_cols=["origin", "malignant_means", "log2_fc", "gene_symbol", "perturbed"],
    pval_col="-log10_pval_adjusted_bh",
    perturbed_col="perturbed",
    facet_col="malignant_means",
    facet_row=None,
    marker_color="perturbed",
)
# 1000 x 500
fig.update_layout(width=1000, height=500)
# png
fig.show(renderer="png", scale=2)
