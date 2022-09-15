import logging

import numpy as np
import pandas as pd
import plotly.express as px
import scipy.stats
from statsmodels.stats.multitest import multipletests

logger = logging.getLogger(__name__)


def compute_stats_for_group(df: pd.DataFrame) -> pd.Series:
    logger.debug("computing stats for group %s", df.name)
    assert len(df.columns) == 2, df.columns
    data_0, data_1 = df.iloc[:, 0], df.iloc[:, 1]
    logger.debug((data_0.shape, data_1.shape))
    sparsity_overall = (df == 0).values.mean()
    results = pd.Series(
        {
            "pval": scipy.stats.mannwhitneyu(data_0, data_1)[1],
            # "pval_ttest": scipy.stats.ttest_ind(data_0, data_1).pvalue,
            "fold_change": data_0.mean() / data_1.mean(),
            "sparsity_overall": sparsity_overall,
        }
    )
    return results


def add_multipletests_stats(df):
    alpha = 0.5  # false discovery rate
    df["-log10_pval"] = -np.log10(df["pval"])
    df["log2_fold_change"] = np.log2(df["fold_change"])
    sign = np.sign(df["log2_fold_change"])
    df["-log10_pval_signed"] = df["-log10_pval"] * sign
    # multiple hypothesis testing with benjamini-hochberg
    df["significant_bh_fdr=0.5"], df["pval_adj_bh"] = multipletests(
        df["pval"], alpha=alpha, method="fdr_bh"
    )[0:2]
    df["-log10_pval_adj_bh"] = -np.log10(df["pval_adj_bh"])
    df["-log10_pval_adj_bh_signed"] = df["-log10_pval_adj_bh"] * sign
    n_signif_results = df["significant_bh_fdr=0.5"].sum()
    df.attrs["pval_threshold_bh"] = (n_signif_results + 1) * alpha / len(df)
    df.attrs["-log10_pval_threshold_bh"] = -np.log10(df.attrs["pval_threshold_bh"])
    return df


def compute_stats(groups):
    df_stats_by_group = groups.apply(compute_stats_for_group)
    df_stats_by_group = df_stats_by_group.reset_index()
    df_stats_by_group = add_multipletests_stats(df_stats_by_group)
    return df_stats_by_group


def make_volcano_figure(df_stats):
    fig = px.scatter(
        df_stats,
        x="log2_fold_change",
        y="-log10_pval",
        color="significant_bh_fdr=0.5",
        hover_name="GeneSymbol",
        hover_data=["pval", "pval_adj_bh", "sparsity_overall"],
    )
    fig.update_layout(
        xaxis_title=r"$\log_{2} [\text{fold change}]$",
        yaxis_title=r"$-\log_{10} [\text{p-value (Mann-Whitney U)}]$",
        legend_title="Significant with FDR=0.5?",
        font=dict(family="Courier New, monospace", color="RebeccaPurple"),
        height=750,
    )
    fig.update_traces(marker=dict(size=3))
    fig.update_layout(legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01))
    fig.add_hline(y=df_stats.attrs["-log10_pval_threshold_bh"])
    return fig
