import logging

import dask.dataframe as dd
import numpy as np
import pandas as pd
import plotly.express as px
import scipy.stats
from google.cloud import bigquery
from statsmodels.stats.multitest import multipletests

from helpers.datasets import load_tcga_skcm_fractions_from_csx

logger = logging.getLogger(__name__)


def make_labels_for_aliquots(df_cell_type_fractions, df_sample_metadata):
    immune_cell_types = ["B", "Macrophage", "NK", "T", "T CD4", "T CD8"]
    df_immune_fraction = (
        df_cell_type_fractions[immune_cell_types]
        .sum(axis="columns")
        .to_frame("immune_fraction")
        .rename_axis(index="aliquot_barcode")
        .reset_index()
        .assign(sample_barcode=lambda row: row["aliquot_barcode"].str[:-12])
    )
    df_sample_metadata = df_sample_metadata[
        ["sample_barcode", "sample_type_name"]
    ].merge(
        df_immune_fraction,
        left_on="sample_barcode",
        right_on="sample_barcode",
        validate="one_to_one",
    )
    immune_quintile = pd.qcut(df_sample_metadata["immune_fraction"], 5, labels=False)
    df_sample_metadata = df_sample_metadata.assign(immune_low=immune_quintile == 0)
    df_sample_metadata = df_sample_metadata.assign(immune_high=immune_quintile == 4)
    return df_sample_metadata


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
    # df["pval_adj_bh_0.05"] = multipletests(df["pval"], method="fdr_bh", alpha=0.05)[1]
    # df["pval_adj_bh_0.10"] = multipletests(df["pval"], method="fdr_bh", alpha=0.1)[1]
    # df["pval_adj_bh_0.20"] = multipletests(df["pval"], method="fdr_bh", alpha=0.2)[1]
    # df["-log10(pval_adj_bh_0.05)"] = -np.log10(df["pval_adj_bh_0.05"])
    # df["-log10(pval_adj_bh_0.10)"] = -np.log10(df["pval_adj_bh_0.10"])
    # df["-log10(pval_adj_bh_0.20)"] = -np.log10(df["pval_adj_bh_0.20"])
    for alpha in [0.05, 0.1, 0.2]:
        field_name = f"pval_adj_fdr={alpha:3.2f}"
        df[field_name] = multipletests(df["pval"], method="fdr_bh", alpha=alpha)[1]
        df[f"-log10_{field_name}"] = -np.log10(df[field_name])
        df[f"-log10_{field_name}_signed"] = df[f"-log10_{field_name}"] * sign
    return df


def get_metastatic_sample_barcodes():
    bqclient = bigquery.Client()
    query_string = """
    SELECT * 
    FROM `isb-cgc-bq.TCGA.biospecimen_gdc_current`
    where project_short_name = "TCGA-SKCM"
        and sample_type_name = "Metastatic"
    order by sample_barcode
    """
    df_tcga_sample_metadata = (
        bqclient.query(query_string).result().to_dataframe(progress_bar_type="tqdm")
    )
    return df_tcga_sample_metadata


def get_tcga_skcm_metastatic_sample_metadata() -> pd.DataFrame:
    df_tcga_skcm_fractions_from_csx = load_tcga_skcm_fractions_from_csx()
    df_tcga_sample_metadata = get_metastatic_sample_barcodes()
    df_sample_metadata = make_labels_for_aliquots(
        df_tcga_skcm_fractions_from_csx, df_tcga_sample_metadata
    )
    return df_sample_metadata


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
        y="-log10_pval_adj_fdr=0.10",
        hover_name="gene_symbol",
        hover_data=["fold_change", "pval"],
    )
    fig.update_xaxes(range=(-12, 12))
    fig.update_yaxes(range=(0, 25))
    fig.update_traces(marker=dict(size=3))
    return fig


def make_scatter_of_signed_pvals(df_gene_stats_1, df_gene_stats_2):
    fig = px.scatter(
        df_gene_stats_1.merge(df_gene_stats_2, on="gene_symbol"),
        x="-log10_pval_adj_fdr=0.10_signed_x",
        y="-log10_pval_adj_fdr=0.10_signed_y",
        trendline="ols",
        hover_name="gene_symbol",
    )
    fig.update_xaxes(range=(-25, 25))
    fig.update_yaxes(range=(-25, 25))
    fig.update_traces(marker=dict(size=3))
    return fig


def make_scatter_of_log2_fold_changes(df_gene_stats_1, df_gene_stats_2):
    fig = px.scatter(
        df_gene_stats_1.merge(df_gene_stats_2, on="gene_symbol"),
        x="log2_fold_change_x",
        y="log2_fold_change_y",
        trendline="ols",
        hover_name="gene_symbol",
    )
    fig.update_xaxes(range=(-12, 12))
    fig.update_yaxes(range=(-12, 12))
    fig.update_traces(marker=dict(size=3))
    return fig


def analyze_gene_significance_overlap(
    df_gene_stats_1: pd.DataFrame,
    df_gene_stats_2: pd.DataFrame,
    top_fraction_cutoff: float,
):
    special_genes_1 = (
        df_gene_stats_1.set_index("gene_symbol")["-log10_pval_adj_fdr=0.10"].rank(
            pct=True
        )
        > 1 - top_fraction_cutoff
    ).rename("significant_in_cohort_1")
    special_genes_2 = (
        df_gene_stats_2.set_index("gene_symbol")["-log10_pval_adj_fdr=0.10"].rank(
            pct=True
        )
        > 1 - top_fraction_cutoff
    ).rename("significant_in_cohort_2")
    special_genes_both = special_genes_1 & special_genes_2
    crosstab = pd.crosstab(special_genes_1, special_genes_2)
    odds_ratio, p_value = scipy.stats.fisher_exact(crosstab)
    return crosstab, odds_ratio, p_value
