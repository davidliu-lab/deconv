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
    log2_fc = -1.0
    AND malignant_means in (
        '0.6,0.8',
        '0.63,0.77',
        '0.65,0.75',
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
# gene-level attributes
sql = f"""
SELECT
    gene_symbol,
    AVG(Malignant) as tpm_avg,
    stddev(Malignant) as tpm_stddev,
    AVG(Malignant = 0) as tpm_zero
FROM
    `deconv.cell_type_geps`
WHERE
    {where_expression}
GROUP BY
    gene_symbol
ORDER BY
    gene_symbol
;
"""
df_gene_attrs = bigquery.Client().query(sql).to_dataframe(progress_bar_type="tqdm")

# %%
# for non-perturbed genes, average p-value by 

# %%
# metrics for different log2_fc values, alpha values
df_metrics = helpers.deg_analysis.classifier_metrics.get_metrics_for_alphas(
    df_gene_stats,
    groupbys=["origin", "malignant_means", "log2_fc", "run_id"],
    score_col="neg_log10_pval_adjusted_bh_signed_directional",
    alphas=[0.1],
).drop(columns=["roc_auc"])

roc_auc = helpers.deg_analysis.classifier_metrics.get_metrics_for_alphas(
    df_gene_stats,
    groupbys=["origin", "malignant_means", "log2_fc", "run_id"],
    score_col="neg_log10_pval_signed_directional",
    alphas=[0.1],
)["roc_auc"]

# %% 
# table
helpers.deg_analysis.displaying_tables.make_score_table_with_stddev(
    df_metrics.query("origin == 'malignant_cibersortx'")["recall"],
    # roc_auc.to_frame("roc_auc").query("origin == 'bulk'")["roc_auc"],
    cmap="Blues",
    index="malignant_means",
    columns="log2_fc",
    vmin=0,
    vmax=2,
    include_n=False,
)

# %%
# improvement
delta = (
    (
        df_metrics
        .query("origin == 'malignant_cibersortx'")
        # drop "origin" index level
        .droplevel("origin")
        ["precision"]
    )
    - (
        df_metrics
        .query("origin == 'bulk'")
        .droplevel("origin")
        ["precision"]
    )
)
helpers.deg_analysis.displaying_tables.make_score_table(
    delta,
    "Blues",
    vmin=0,
    vmax=.3,
)

# %%
# ROC AUC
helpers.deg_analysis.displaying_tables.make_score_table_with_stddev(
    df_metrics.query("origin == 'malignant_cibersortx'")["recall"],
    cmap="Blues",
    index="malignant_means",
    columns="log2_fc",
    vmin=0,
    vmax=1,
    include_n=False,
)


# %%
# volcano plots for specific cell

df_curves = helpers.deg_analysis.classifier_metrics.get_curves_with_all_pvals(
    df_gene_stats.query((
        "malignant_means == '0.77,0.63'"
        " AND log2_fc == 1.0"
        " and run_id == 0"
    )),
    groupby_cols=["origin", "malignant_means", "log2_fc", "run_id"],
    score_col="neg_log10_pval_signed_directional",
)


# %%
# ROC curves comparing
