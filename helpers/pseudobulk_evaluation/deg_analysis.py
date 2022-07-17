import logging
import pandas as pd
import scipy.stats
from statsmodels.stats.multitest import multipletests
import numpy as np


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


def add_immune_groups_to_rnaseq(ddf_bulk_rnaseq, df_sample_metadata):
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
    logger.debug(df.dtypes)
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
    sign = df["log2_fold_change"].map(lambda x: 1 if x > 0 else -1)
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
