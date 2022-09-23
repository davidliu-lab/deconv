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
    data_baseline, data_other = df.iloc[:, 0], df.iloc[:, 1]
    logger.debug((data_baseline.shape, data_other.shape))
    sparsity_overall = (df == 0).values.mean()
    results = pd.Series(
        {
            "pval": scipy.stats.mannwhitneyu(data_baseline, data_other)[1],
            # "pval_ttest": scipy.stats.ttest_ind(data_0, data_1).pvalue,
            "fold_change": data_other.mean() / data_baseline.mean(),
            "sparsity_overall": sparsity_overall,
        }
    )
    return results


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
        df.attrs[f"-log10_{pval_threshold_str}"] = -np.log10(df.attrs["pval_threshold_bh"])
    return df


def compute_stats(groups) -> pd.DataFrame:
    df_stats_by_group = groups.apply(compute_stats_for_group)
    df_stats_by_group = df_stats_by_group.reset_index()
    df_stats_by_group = add_multipletests_stats(df_stats_by_group)
    return df_stats_by_group


def make_volcano_figure(df_stats: pd.DataFrame):
    fig = px.scatter(
        df_stats,
        x="log2_fold_change",
        y="-log10_pval",
        color="perturbed",
        hover_name="gene_symbol",
        hover_data=["pval", "sparsity_overall", "perturbation_factor"],
    )
    fig.update_layout(
        xaxis_title=r"$\log_{2} [\text{fold change}]$",
        yaxis_title=r"$-\log_{10} [\text{p-value (Mann-Whitney U)}]$",
        legend_title="Perturbed?",
        font=dict(family="Courier New, monospace", color="RebeccaPurple"),
        height=750,
    )
    fig.update_traces(marker=dict(size=3))
    fig.update_layout(legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01))
    fig.add_hline(y=df_stats.attrs["-log10_pval_threshold_fdr=0.10"])
    fig.add_hline(y=df_stats.attrs["-log10_pval_threshold_fdr=0.25"])
    return fig
