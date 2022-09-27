import logging
import warnings

import numpy as np
import pandas as pd
from pandas.core.groupby import SeriesGroupBy
import plotly.express as px
import scipy.stats
from statsmodels.stats.multitest import multipletests

logger = logging.getLogger(__name__)


def add_multipletests_stats(df: pd.DataFrame) -> pd.DataFrame:
    alphas = [0.1, 0.25]  # false discovery rates
    df["-log10_pval"] = -np.log10(df["pval"])
    df["log2_fold_change"] = np.log2(df["fold_change"])
    sign = np.sign(df["log2_fold_change"])
    df["-log10_pval_signed"] = df["-log10_pval"] * sign
    # multiple hypothesis testing with benjamini-hochberg
    for alpha in alphas:
        significance_column = f"significant_bh_fdr={alpha:.2f}"
        df[significance_column] = multipletests(
            df["pval"], alpha=alpha, method="fdr_bh"
        )[0]
        n_signif_results = df[significance_column].sum()
        pval_threshold_str = f"pval_threshold_fdr={alpha:.2f}"
        df.attrs[pval_threshold_str] = (n_signif_results + 1) * alpha / len(df)
        df.attrs[f"-log10_{pval_threshold_str}"] = -np.log10(
            df.attrs[pval_threshold_str]
        )
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
            }
        )
        return results

    df_stats_by_group = series_groupby.apply(compute_stats_for_group)
    df_stats_by_group = df_stats_by_group.unstack(-1)
    df_stats_by_group = df_stats_by_group.reset_index()
    df_stats_by_group = add_multipletests_stats(df_stats_by_group)
    return df_stats_by_group


def make_volcano_figure(df_stats: pd.DataFrame, perturbed: bool = False):
    fig = px.scatter(
        df_stats,
        x="log2_fold_change",
        y="-log10_pval",
        color="perturbed" if perturbed else None,
        hover_name="gene_symbol",
        hover_data=["pval", "sparsity_overall"],
    )
    fig.update_layout(
        xaxis_title=r"$\log_{2} [\text{fold change}]$",
        yaxis_title=r"$-\log_{10} [\text{p-value (Mann-Whitney U)}]$",
        legend_title="Perturbed?" if perturbed else None,
        font=dict(family="Courier New, monospace", color="RebeccaPurple"),
        height=750,
    )
    fig.update_traces(marker=dict(size=3))
    fig.update_layout(legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01))
    fig.add_hline(y=df_stats.attrs["-log10_pval_threshold_fdr=0.10"])
    fig.add_hline(y=df_stats.attrs["-log10_pval_threshold_fdr=0.25"])
    return fig
