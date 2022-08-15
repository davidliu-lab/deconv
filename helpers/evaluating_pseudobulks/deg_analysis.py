import logging

import dask.dataframe as dd
import numpy as np
import pandas as pd
import plotly.express as px
import scipy.stats
from statsmodels.stats.multitest import multipletests

logger = logging.getLogger(__name__)


def compute_gene_stats_with_immune_groups(ddf_bulk_rnaseq, df_sample_metadata):
    df = ddf_bulk_rnaseq.join(
        df_sample_metadata.set_index("aliquot_barcode")[
            ["immune_fraction", "immune_low", "immune_high"]
        ],
        on="aliquot_barcode",
    ).compute()
    df_gene_data = df.groupby("gene_symbol").apply(compute_stats_single_gene)
    df_gene_data = df_gene_data.reset_index()
    return df_gene_data


def compute_stats_single_gene(df: pd.DataFrame) -> pd.Series:
    """
    Compute differential expression stats for a single gene

    :param df: pandas.DataFrame, rows are samples
    """
    immune_low_tpm = df[df["immune_low"]]["tpm"]
    immune_high_tpm = df[df["immune_high"]]["tpm"]
    pval = scipy.stats.mannwhitneyu(immune_high_tpm, immune_low_tpm)[1]
    fold_change = immune_high_tpm.mean() / immune_low_tpm.mean()
    return pd.Series({"pval": pval, "fold_change": fold_change})


def process_gene_level_results(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add adjusted p-values for multiple testing using Benjamini-Hochberg procedure, and other fields

    :param df: pandas.DataFrame, rows are genes
    """
    df = df.copy()
    df["-log10_pval"] = -np.log10(df["pval"])
    df["log2_fold_change"] = np.log2(df["fold_change"])
    sign = np.sign(df["log2_fold_change"])
    df["-log10_pval_signed"] = df["-log10_pval"] * sign

    # benjamini-hochberg adjustments
    ## using scipy.stats.multipletests
    df["pval_adj_bh"] = multipletests(df["pval"], method="fdr_bh")[1]
    df["-log10_pval_adj_bh"] = -np.log10(df["pval_adj_bh"])
    df["-log10_pval_adj_bh_signed"] = df["-log10_pval_adj_bh"] * sign

    ## manually using pandas.DataFrame.rank
    df["pval_rank_pandas"] = df["pval"].rank(method="min")
    df["pval_adj_rank_pandas"] = df["pval"] * len(df) / df["pval_rank_pandas"]
    df["-log10_pval_adj_rank_pandas"] = -np.log10(df["pval_adj_rank_pandas"])

    for alpha in [0.05, 0.1, 0.2]:
        field_name = f"significant_bh_fdr={alpha:3.2f}"
        df[field_name] = multipletests(df["pval"], method="fdr_bh", alpha=alpha)[0]
        logger.debug("computing B-H significance test manually for alpha=%3.2f", alpha)
        df[f"{field_name}_rank_pandas"] = False
        for row_index in df["pval"].sort_values().index:
            if df.loc[row_index, "pval_adj_rank_pandas"] <= alpha:
                df.loc[row_index, f"{field_name}_rank_pandas"] = True
            else:
                logger.debug(
                    "stopping manual B-H search at rank %d",
                    df.loc[row_index, "pval_rank_pandas"],
                )
                break
    return df


def compute_all_deg_results(
    ddf_bulk_rnaseq: dd.DataFrame, df_sample_metadata: pd.DataFrame
) -> pd.DataFrame:
    logger.debug("computing stats")
    df_gene_stats = compute_gene_stats_with_immune_groups(
        ddf_bulk_rnaseq, df_sample_metadata
    )
    logger.debug("creating derived columns (e.g. -log10_pval...)")
    df_gene_stats = process_gene_level_results(df_gene_stats)
    return df_gene_stats


def make_volcano_plot(df_gene_stats):
    fig = px.scatter(
        df_gene_stats,
        x="log2_fold_change",
        y="-log10_pval_adj_bh",
        color="significant_bh_fdr=0.10",
        color_discrete_map={True: "orangered", False: "royalblue"},
        hover_name="gene_symbol",
        hover_data=["fold_change", "pval"],
    )
    fig.update_xaxes(range=(-12, 12))
    fig.update_yaxes(range=(0, 25))
    fig.update_traces(marker=dict(size=2.5))
    fig.update_layout(
        legend=dict(
            yanchor="top", y=0.99, xanchor="left", x=0.01, itemsizing="constant"
        ),
    )
    return fig


def make_scatter_of_signed_pvals(df_gene_stats_merged):
    df = df_gene_stats_merged
    df["significant_in"] = df.apply(signif_in, axis=1)
    fig = px.scatter(
        df,
        x="-log10_pval_adj_bh_signed_x",
        y="-log10_pval_adj_bh_signed_y",
        color="significant_in",
        # color_discrete_map={True: "orangered", False: "royalblue"},
        color_discrete_sequence=px.colors.qualitative.G10,
        hover_name="gene_symbol",
    )
    fig.update_xaxes(range=(-25, 25))
    fig.update_yaxes(range=(-25, 25))
    fig.update_traces(marker=dict(size=2.5))
    fig.update_layout(
        legend=dict(
            yanchor="top", y=0.99, xanchor="left", x=0.01, itemsizing="constant"
        ),
    )
    return fig


def signif_in(row):
    x = row["significant_bh_fdr=0.10_x"]
    y = row["significant_bh_fdr=0.10_y"]
    if x & y:
        return "both"
    if x:
        return "x only"
    if y:
        return "y only"
    return "neither"


def make_scatter_of_log2_fold_changes(df_gene_stats_1, df_gene_stats_2):
    df = df_gene_stats_1.merge(df_gene_stats_2, on="gene_symbol")

    df["significant_in"] = df.apply(signif_in, axis=1)
    fig = px.scatter(
        df,
        x="log2_fold_change_x",
        y="log2_fold_change_y",
        symbol="significant_in",
        # trendline="ols",
        color="significant_in",
        color_discrete_map={
            "both": "orangered",
            "x only": "orchid",
            "y only": "yellow",
            "neither": "royalblue",
        },
        symbol_map={
            "both": "square",
            "x only": "circle",
            "y only": "circle",
            "neither": "circle",
        },
        hover_name="gene_symbol",
    )
    fig.update_xaxes(range=(-10, 10))
    fig.update_yaxes(range=(-10, 10))
    fig.update_traces(marker=dict(size=2.5))
    fig.update_layout(
        legend=dict(
            yanchor="top", y=0.99, xanchor="left", x=0.01, itemsizing="constant"
        ),
    )
    return fig


def analyze_gene_significance_overlap(
    special_genes_1,
    special_genes_2,
):
    crosstab = pd.crosstab(special_genes_1, special_genes_2)
    odds_ratio, p_value = scipy.stats.fisher_exact(crosstab)
    return crosstab, odds_ratio, p_value
