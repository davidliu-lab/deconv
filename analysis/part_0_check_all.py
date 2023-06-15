# %%[markdown]
# # Analysis - check all results
# goals
# - sanity check results

# %%
# imports, logging, configuration
import logging

import duckdb

import helpers
from helpers.deg_analysis import displaying_tables, plotting_curves, plotting_volcanos
from helpers.deg_analysis.classifier_metrics import calculate_all_curves
from helpers.deg_analysis.postprocessing_gene_stats_fields import add_more_pval_fields
from helpers.running_cibersortx.loading_results import (
    get_arrow_dataset_for_deg_analysis_results,
)

logging.basicConfig(format=helpers.logging.FORMAT, level="INFO")
logger = logging.getLogger(__name__)
logger.setLevel("DEBUG")
logging.getLogger("helpers.deg_analysis").setLevel("DEBUG")

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
    origin = 'malignant_cibersortx'
    --AND malignant_means in ('0.6,0.8', 'None,None', '0.8,0.6')
    --AND log2_fc in (-1.0, 0.0, 1.0)
    --AND run_id in (0, 1)
    --AND gene_symbol like '[A-K]%'  -- start with A-K
;
"""
df_gene_stats = duckdb.sql(query_text).df()  # query with "duckdb:///:default:"
df_gene_stats = add_more_pval_fields(df_gene_stats)


# %%
# compute curves

df_curves_pval_signed_directional = calculate_all_curves(
    df_gene_stats,
    score_col="-log10_pval_signed_directional",
    perturbed_col="perturbed",
)
df_curves_pval_adjusted_bh_signed_directional = calculate_all_curves(
    df_gene_stats,
    score_col="-log10_pval_adjusted_bh_signed_directional",
    perturbed_col="perturbed",
)
_, _, df_scores = helpers.deg_analysis.classifier_metrics.compute_all_curves_and_metrics(
    df_gene_stats, signed_directional=True
)


# %%
displaying_tables.make_score_table_with_stddev(df_scores["roc_auc"], "RdBu")


# %%
displaying_tables.make_score_table_with_stddev(df_scores["precision"], "RdBu")

# %%
fig = plotting_curves.plot_roc_curves(df_curves_pval_signed_directional.reset_index())
fig.update_layout(width=1000, height=1000)
fig.show(renderer="png", scale=2)

# %% [markdown]
# check...
# - are all values of "-log10_pval_signed_directional" in the dataframe df_gene_stats?
# - can i do a join to get corresponding values of "-log10_pval_adjusted_bh_signed_directional"?

logger.setLevel("DEBUG")

index_names = df_curves_pval_signed_directional.index.names
# logger.debug("index_names: %s", index_names)
# join_keys = list(index_names)
# logger.debug("join_keys: %s", join_keys)
# cols_from_df_gene_stats = join_keys + ["-log10_pval_adjusted_bh_signed_directional"]
# logger.debug("cols_from_df_gene_stats: %s", cols_from_df_gene_stats)
new_col = "-log10_pval_adjusted_bh_signed_directional"
other = df_gene_stats.groupby(index_names)[new_col].first()

df_curves_with_both_pvals = df_curves_pval_signed_directional.join(
    other=other,
    how="left",
    validate="one_to_one",
)

print(df_curves_pval_signed_directional.shape)
print(df_curves_with_both_pvals.shape)

# %% [markdown]
# # table of ROC AUC values
displaying_tables.make_score_table_with_stddev(df_scores["roc_auc"], "RdBu")


# %% [markdown]
# # table of ROC AUC values
displaying_tables.make_score_table_with_stddev(df_scores["precision"], "RdBu")


# %% [markdown]
# # plots of ROC curves
fig = plotting_curves.plot_roc_curves(df_curves_with_both_pvals.reset_index())
fig.update_layout(width=1000, height=1000)
fig.show(renderer="png", scale=2)

# %% [markdown]
# # plots of precision-recall curves
fig = plotting_curves.plot_precision_recall_curves(df_curves_with_both_pvals.reset_index())
fig.update_layout(width=1000, height=1000)
fig.show(renderer="png", scale=2)


# %% [markdown]
# ## volcano plots
fig = plotting_volcanos.make_volcano_grid_scatter(
    df_gene_stats,
    perturbed_col="perturbed",
    aggregate_by=["origin", "malignant_means", "log2_fc", "gene_symbol", "perturbed"],
)
# 1000 x 1500
fig.update_layout(width=1000, height=1500)
