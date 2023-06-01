import logging
import warnings

import numpy as np
import pandas as pd
import scipy.stats
from pandas.core.groupby.generic import SeriesGroupBy
from statsmodels.stats.multitest import multipletests

logger = logging.getLogger(__name__)


def calculate_pval_threshold(test_was_fdr_significant: pd.Series, alpha: float) -> float:
    n_signif_results = test_was_fdr_significant.sum()
    pval_threshold = (n_signif_results + 1) * alpha / len(test_was_fdr_significant)
    return pval_threshold


def add_multipletests_stats(df: pd.DataFrame) -> pd.DataFrame:
    alphas = [0.1, 0.25]  # false discovery rates
    df["-log10_pval"] = -np.log10(df["pval"])
    df["log2_fold_change"] = np.log2(df["fold_change"])
    sign = np.sign(df["log2_fold_change"])
    df["-log10_pval_signed"] = df["-log10_pval"] * sign
    # multiple hypothesis testing with benjamini-hochberg
    for alpha in alphas:
        significance_column = f"significant_bh_fdr={alpha:.2f}"
        reject, pvals_corrected = multipletests(df["pval"], alpha=alpha, method="fdr_bh")[0:2]
        df[significance_column] = reject
        df["pval_adjusted_bh"] = pvals_corrected
        n_signif_results = df[significance_column].sum()
        pval_threshold_str = f"pval_threshold_fdr={alpha:.2f}"
        df.attrs[pval_threshold_str] = (n_signif_results + 1) * alpha / len(df)
        df.attrs[f"-log10_{pval_threshold_str}"] = -np.log10(df.attrs[pval_threshold_str])
    return df


def compute_stats(
    series_groupby: SeriesGroupBy,
    group_col: str,
    group_1: str,
    group_2: str,
) -> pd.DataFrame:
    def compute_stats_for_group(series: pd.Series) -> pd.Series:
        logger.debug("computing stats for group %s", series.name)
        series_groupby = series.groupby(group_col)
        series_1 = series_groupby.get_group(group_1)
        series_2 = series_groupby.get_group(group_2)
        logger.debug("shapes: %s", (series_1.shape, series_2.shape))
        sparsity_overall = np.mean(series == 0)
        logger.debug("sparsity overall: %s", sparsity_overall)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # ignore divide by zero warnings
            fold_change = series_2.mean() / series_1.mean()
        results = pd.Series(
            {
                "pval": scipy.stats.mannwhitneyu(series_1, series_2)[1],
                # "pval_ttest": scipy.stats.ttest_ind(series_1, series_2).pvalue,
                "fold_change": fold_change,
                "sparsity_overall": sparsity_overall,
                "expression_mean_1": series_1.mean(),
                "expression_mean_2": series_2.mean(),
            }
        )
        return results

    df_stats_by_group = series_groupby.apply(compute_stats_for_group)
    df_stats_by_group = df_stats_by_group.unstack(-1)
    df_stats_by_group = df_stats_by_group.reset_index()
    df_stats_by_group = add_multipletests_stats(df_stats_by_group)
    return df_stats_by_group
