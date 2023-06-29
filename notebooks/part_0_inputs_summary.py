# %%
# imports, logging, configuration
import logging

import duckdb
import numpy as np
import pandas as pd
import plotly.express as px
from IPython.core.interactiveshell import InteractiveShell
from google.cloud import bigquery

import helpers
from helpers.deg_analysis import plotting_curves, plotting_utils
from helpers.deg_analysis.loading_results import (
    get_arrow_dataset_for_deg_analysis_results,
)
from helpers.deg_analysis.plotting_curves import plot_curves, plot_metric_by_threshold
from helpers.deg_analysis.plotting_utils import add_fdr_lines
from helpers.deg_analysis.plotting_volcanos import make_volcano_grid_scatter
from helpers.deg_analysis.postprocessing_gene_stats_fields import add_more_pval_fields
from helpers.logging.configuring_logging import FORMAT, EasternTimeFormatter

logger = logging.getLogger(__name__)
logger.setLevel("DEBUG")
logging.getLogger("helpers").setLevel("INFO")
logging.getLogger("google.cloud").setLevel("INFO")
handler = logging.StreamHandler()
handler.setFormatter(EasternTimeFormatter(FORMAT, datefmt="%I:%M:%S %p"))
logging.getLogger().handlers = [handler]

# %%
logger.debug("test")

# %%
# ipython magics
if InteractiveShell.initialized():
    ipython = InteractiveShell.instance()
    ipython.run_line_magic("load_ext", "autoreload")
    ipython.run_line_magic("autoreload", "2")
    ipython.run_line_magic("sql", "duckdb:///:default:")

# %%
where_expression = """
    --malignant_means = '0.85,0.55'
    malignant_means = '0.55,0.85'
    --malignant_means = '0.65,0.75'
    --malignant_means = 'None,None'  -- no perturbation
    AND log2_fc in (1.0)
    AND run_id = 0
"""

# %%
# load gene stats

bq_client = bigquery.Client()

sql = f"""
SELECT *
FROM `deconv.deg_analysis_results`
WHERE {where_expression}
;
"""
df_gene_stats = bq_client.query(sql).to_dataframe(progress_bar_type="tqdm")

# %%
len(df_gene_stats)

# %%
fractions = helpers.running_cibersortx.reading_input_files.get_arrow_dataset(
    "gs://liulab/differential_composition_and_expression/20230616_03h34m20s/fractions/"
)
df_fractions = duckdb.sql(f"""
SELECT *
FROM fractions
WHERE {where_expression}
""").to_df()
df_fractions

# %%
# stacked bar of fractions
fig = px.bar(
    df_fractions.sort_values(["name", "Malignant"]),
    x="sample_id",
    y=list(helpers.cell_type_naming.nice_to_weirds.keys()),
    title="Cell type proportions, cohorts A and B"
)
# hide x-axis labels
fig.update_xaxes(showticklabels=False)
# set title: "Cell type proportions, cohorts A and B"
fig.show(renderer="png", scale=2, width=800, height=500)

# %%
# volcano plot: DGE for bulk RNA-seq

fig = px.scatter(
    df_gene_stats.query("origin == 'bulk'"),
    x="log2_fold_change",
    y="neg_log10_pval_adjusted_bh",
    color="perturbed",
    # facet_col="origin",
    hover_name="gene_symbol",
)
fig = helpers.deg_analysis.plotting_volcanos.format_volcano_figure(fig)
# fig.update_xaxes(range=[-4, 4])
# fig.update_layout(title="DEG analysis of bulk RNA-seq and CIBERSORTx malignant-specific outputs")
# set title: "Differential expression: bulk RNA-seq"
fig.update_layout(title="Differential expression: bulk RNA-seq")
fig.show(renderer="png", scale=2, width=500, height=500)


# %%
cell_type_geps = helpers.running_cibersortx.reading_input_files.get_arrow_dataset(
    "gs://liulab/differential_composition_and_expression/20230616_03h34m20s/cell_type_geps/"
)
df_cell_type_geps = duckdb.sql(f"""
SELECT *
FROM cell_type_geps
WHERE {where_expression}
""").to_df()

df_cell_type_geps

# %%
# compute stats for malignant cell type GEPs
logging.getLogger("helpers.deg_analysis").setLevel("INFO")

def get_malignant_stats(df_geps):
    groupby_cols = ["malignant_means", "log2_fc", "run_id", "gene_symbol"]
    df_geps = df_geps.set_index(groupby_cols + ["name"])
    series_groupby = df_geps.groupby(groupby_cols)["Malignant"]
    df = helpers.deg_analysis.stats_testing_with_fdr.compute_stats(
        series_groupby,
        group_col="name",
        group_1="a",
        group_2="b",
    )
    df.rename(columns=lambda x: x.replace("-log10", "neg_log10"), inplace=True)
    df["neg_log10_pval_adjusted_bh"] = df["pval_adjusted_bh"].apply(lambda x: -np.log10(x))
    genes_perturbed = pd.read_csv("gs://liulab/differential_composition_and_expression/20230616_03h34m20s/genes_perturbed.csv")
    df["perturbed"] = (df["log2_fc"] != 0.0) & df["gene_symbol"].apply(lambda x: x in genes_perturbed["gene_symbol"].values)
    return df

malignant_stats = get_malignant_stats(df_cell_type_geps)

# %%

# malignant_stats = malignant_stats.assign(**{
#     "neg_log10_pval_adjusted_bh": lambda x: -np.log10(x["pval_adjusted_bh"])
# })
assert "neg_log10_pval_adjusted_bh" in malignant_stats.columns, list(malignant_stats.columns)

# %%
# volcano plot: DGE for malignant cell type GEP
fig = px.scatter(
    malignant_stats,
    x="log2_fold_change",
    y="neg_log10_pval_adjusted_bh",
    color="perturbed",
    hover_name="gene_symbol",
)
fig = helpers.deg_analysis.plotting_volcanos.format_volcano_figure(fig)
# set title: "Differential expression: malignant cell type GEPs"
fig.update_layout(title="Differential expression: true malignant-specific GEP")
fig.show(renderer="png", scale=2, width=500, height=500)


# %%
# genes where, in malignant, CIBERSORTx q-val == 0.0
genes_cibersortx_qval_0 = (
    df_gene_stats[
        (df_gene_stats["neg_log10_pval_adjusted_bh"] == 0.0)
        & (df_gene_stats["origin"] == "malignant_cibersortx")
    ]["gene_symbol"]
    .tolist()
)
print(len(genes_cibersortx_qval_0), "genes")

# %%
# compare log2_fold_change in bulk, comparing gene in genes_cibersortx_qval_0 vs not


df_to_plot = (
    df_gene_stats
    [df_gene_stats["origin"] == "bulk"]
    [["gene_symbol", "log2_fold_change"]]
    .assign(
        in_cibersortx_qval_0=lambda df: df["gene_symbol"].isin(genes_cibersortx_qval_0)
    )
)
df_to_plot
px.box(
    df_to_plot,
    x="in_cibersortx_qval_0",
    y="log2_fold_change",
    hover_name="gene_symbol",
    width=500,
    height=500,
)

# %%
# find genes where log2_fold_change is neither nan nor inf
genes_nan = (
    df_gene_stats
    [df_gene_stats["log2_fold_change"].isna()]
    ["gene_symbol"]
    .unique()
    .tolist()
)
genes_inf = (
    df_gene_stats
    [df_gene_stats["log2_fold_change"].isin([float("inf"), float("-inf")])]
    ["gene_symbol"]
    .unique()
    .tolist()
)
print(len(genes_nan), "genes with nan log2_fold_change")
print(len(genes_inf), "genes with inf log2_fold_change")

# %%
# filter out inf or nan
df_gene_stats = df_gene_stats[
    ~df_gene_stats["gene_symbol"].isin(genes_nan + genes_inf)
]
len(df_gene_stats)

# %%
# scatter of p-values, bulk vs CIBERSORTx
gene_level_attributes = [
    "neg_log10_pval_adjusted_bh",
    "log2_fold_change",
]
df_ = (
    df_gene_stats
    [["gene_symbol", "origin", "perturbed"] + gene_level_attributes]
    .set_index(["gene_symbol", "perturbed", "origin"])
    .unstack("origin")["neg_log10_pval_adjusted_bh"]
    .reset_index()
)
fig = px.scatter(
    df_,
    x="bulk",
    y="malignant_cibersortx",
    facet_col="perturbed",
    color="perturbed",
    hover_name="gene_symbol",
    # trendline="ols",
)
import plotly.graph_objects as go
# add line y=x
fig.add_trace(
    go.Scatter(
        x=[0, 10],
        y=[0, 10],
        mode="lines",
        line=go.scatter.Line(color="gray", dash="dash"),
        showlegend=False,
    )
)
fig.update_layout(
    title="P-values of genes, comparing bulk and CIBERSORTx-inferred malignant GEP",
    width=1000,
    height=500,
)
# %%
# distributions of log2_fold_change in bulk, comparing CIBERSORTx q-val == 0.0 vs != 0.0
gene_level_attributes = [
    "neg_log10_pval_adjusted_bh",
    "log2_fold_change",
]
df_ = (
    df_gene_stats
    [["gene_symbol", "origin", "perturbed"] + gene_level_attributes]
    .set_index(["gene_symbol", "perturbed", "origin"])
    .unstack("origin")
    .reset_index()
)
# show hist of log2_fold_change
df_to_plot = df_.copy()
df_["log2_fold_change_in_bulk"] = df_[("log2_fold_change", "bulk")]
df_["CIBERSORTx q-val == 0.0"] = df_[("neg_log10_pval_adjusted_bh", "malignant_cibersortx")] == 0.0
fig = px.histogram(
    df_,
    x="log2_fold_change_in_bulk",
    color="CIBERSORTx q-val == 0.0",
    width=1000,
    height=500,
)
fig = px.ecdf(
    df_,
    x="log2_fold_change_in_bulk",
    color="CIBERSORTx q-val == 0.0",
    width=1000,
    height=500,
)
fig


# %%
df_["log2_fold_change_in_bulk"].isna().sum()
# %%
# mean and stdev of log2_fold_change in bulk, comparing CIBERSORTx q-val == 0.0 vs != 0.0
(
    df_
    [["CIBERSORTx q-val == 0.0", "log2_fold_change_in_bulk"]]
    .groupby("CIBERSORTx q-val == 0.0")
    ["log2_fold_change_in_bulk"]
    .agg(["count", "mean", "std"])
)


# %%
# average log2_fold_change in bulk, comparing when neg_log10_pval_adjusted_bh is 0.0 vs > 0.0
# drop nans
df_ = df_.dropna(subset=["log2_fold_change_in_bulk"])
df_.groupby("CIBERSORTx q-val == 0.0")["log2_fold_change_in_bulk"].mean()

# %%
df_

# %%
# what gene-level features predict genes that have low neg_log10_pval_adjusted_bh in ciberortx, but high in bulk?

# %%
# compute curves
df_curves = helpers.deg_analysis.classifier_metrics.get_metrics_for_alphas(
    df_gene_stats,
    groupbys=["origin", "malignant_means", "run_id"],
    score_col="neg_log10_pval_adjusted_bh_signed_directional",
)
df_curves[["fpr", "precision", "recall"]].rename_axis(columns="metrics").unstack("origin").stack("metrics")


# %%
# histogram of log2_fold_change for each origin
fig = px.histogram(
    df_gene_stats,
    x="log2_fold_change",
    facet_col="origin",
    nbins=500,
    histnorm="percent",
)
fig.update_layout(title="Distribution of between-group fold changes in gene expression")
fig.update_xaxes(range=[-10, 10])
fig.update_yaxes(range=[0, 100])
fig.update_layout(width=800, height=400)
fig.show(renderer="png", scale=2)
