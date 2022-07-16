import pandas as pd
import scipy.stats
import numpy as np


def compute_stats_all_genes_with_fdr_adjustments(df):
    """
    Compute DEG stats for all genes with FDR adjusted p-values

    :param df: pandas.DataFrame, rows are samples and genes
    """
    immune_low, immune_high = df[df["immune_low"]], df[df["immune_high"]]
    
    pval = scipy.stats.mannwhitneyu(
        immune_high["tpm"].values, immune_low["tpm"].values
    )[1]
    fold_change = immune_high["tpm"].mean() / immune_low["tpm"].mean()
    return pd.Series(
        dict(
            pval=pval,
            fold_change=fold_change,
        )
    )


def compute_stats_single_gene(df):
    """
    Compute differential expression stats for a single gene

    :param df: pandas.DataFrame, rows are samples
    """
    immune_low, immune_high = df[df["immune_low"]], df[df["immune_high"]]
    pval = scipy.stats.mannwhitneyu(
        immune_high["tpm"].values, immune_low["tpm"].values
    )[1]
    fold_change = immune_high["tpm"].mean() / immune_low["tpm"].mean()
    return pd.Series(
        dict(
            pval=pval,
            fold_change=fold_change,
        )
    )


def compute_stats(df):
    """
    Compute some differential expression statistics for evaluation pseudobulk RNA-seq samples
    """
    immune_low = df[df["immune_low"]]
    immune_high = df[df["immune_high"]]
    pval_pseudo = scipy.stats.mannwhitneyu(
        immune_high["tpm_pseudo"].values, immune_low["tpm_pseudo"].values
    )[1]
    foldchange_pseudo = (
        immune_high["tpm_pseudo"].mean() / immune_low["tpm_pseudo"].mean()
    )
    neglog10pval_pseudo = -np.log10(pval_pseudo)
    log2foldchange_pseudo = np.log2(foldchange_pseudo)

    pval_real = scipy.stats.mannwhitneyu(
        immune_high["tpm_tcga_skcm"].values, immune_low["tpm_tcga_skcm"].values
    )[1]
    foldchange_real = (
        immune_high["tpm_tcga_skcm"].mean() / immune_low["tpm_tcga_skcm"].mean()
    )
    neglog10pval_real = -np.log10(pval_real)
    log2foldchange_real = np.log2(foldchange_real)

    return pd.Series(
        dict(
            pval_pseudo=pval_pseudo,
            foldchange_pseudo=foldchange_pseudo,
            log2foldchange_pseudo=log2foldchange_pseudo,
            neglog10pval_pseudo=neglog10pval_pseudo,
            signedneglog10pval_pseudo=(
                neglog10pval_pseudo * np.sign(log2foldchange_pseudo)
            ),
            pval_real=pval_real,
            foldchange_real=foldchange_real,
            log2foldchange_real=log2foldchange_real,
            neglog10pval_real=neglog10pval_real,
            signedneglog10pval_real=(neglog10pval_real * np.sign(log2foldchange_real)),
        )
    )
