import logging
import pandas as pd
import scipy.stats
from statsmodels.stats.multitest import multipletests
import numpy as np


logger = logging.getLogger(__name__)


def compute_stats_single_gene(df: pd.DataFrame) -> pd.Series:
    """
    Compute differential expression stats for a single gene

    :param df: pandas.DataFrame, rows are samples
    """
    logger.debug(df.info())
    logger.debug(df.dtypes)
    logger.debug(repr(df))
    immune_low_tpm = df[df["immune_low"]]["tpm"]
    immune_high_tpm = df[df["immune_high"]]["tpm"]
    pval = scipy.stats.mannwhitneyu(immune_high_tpm, immune_low_tpm)[1]
    fold_change = immune_high_tpm.mean() / immune_low_tpm.mean()
    return pd.Series({"pval": pval, "fold_change": fold_change})


def adjust_pvals(df: pd.DataFrame) -> pd.DataFrame:
    """
    Adjust p-values for multiple testing using Benjamini-Hochberg procedure

    :param df: pandas.DataFrame, rows are genes
    """
    df["pval_adj_bh_0.05"] = multipletests(df["pval"], method="fdr_bh", alpha=0.05)[1]
    df["pval_adj_bh_0.10"] = multipletests(df["pval"], method="fdr_bh", alpha=0.1)[1]
    df["pval_adj_bh_0.20"] = multipletests(df["pval"], method="fdr_bh", alpha=0.2)[1]
    return df
