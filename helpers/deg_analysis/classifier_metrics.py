"""
Functions for computing classifier metrics (ROC, PR, etc.) from a dataframe of
gene-level statistics.
"""
import logging
import warnings

import numpy as np
import pandas as pd
from pandas.core.groupby.generic import DataFrameGroupBy
from sklearn.metrics import roc_curve

logger = logging.getLogger(__name__)


def get_metrics_for_threshold(
    gene_stats_groupby: DataFrameGroupBy,
    threshold: float = -np.log10(0.1),
    score_col: str = "-log10_pval_adjusted_bh_signed_directional",
):
    """Get TPR, FPR, precision, recall, etc. for a given threshold and score column."""

    def compute_stuff(df) -> pd.Series:
        y_true = df["perturbed"]
        y_pred = df[score_col] >= threshold
        data = {
            "n": len(df),
            "n_true": y_true.sum(),
            "n_false": len(df) - y_true.sum(),
            "fp": y_pred.sum(),
            "fpr": y_pred.sum() / len(df),
        }
        try:
            assert y_true.sum() != 0
            data["tp"] = (y_true & y_pred).sum()
            data["tpr"] = (y_true & y_pred).sum() / y_true.sum()
            data["precision"] = (y_true & y_pred).sum() / y_pred.sum()
            data["recall"] = (y_true & y_pred).sum() / y_true.sum()
        # except assertion error or zero division error
        except (AssertionError, ZeroDivisionError):
            data["tp"] = np.nan
            data["tpr"] = np.nan
            data["precision"] = np.nan
            data["recall"] = np.nan
        return pd.Series(data)

    return gene_stats_groupby.apply(compute_stuff)


def get_metrics_for_alphas(
    df_gene_stats,
    groupbys: list[str] = ["malignant_means", "run_id"],
    alphas: list[float] = [0.05, 0.1, 0.25],
):
    return (
        pd.concat(
            {
                alpha: get_metrics_for_threshold(
                    df_gene_stats.groupby(groupbys),
                    threshold=-1.0 * np.log10(alpha),
                    score_col="-log10_pval_adjusted_bh_signed_directional",
                )
                for alpha in alphas
            },
            names=["alpha"],
        )
        .reorder_levels(groupbys + ["alpha"])
        .sort_index()
    )


def get_curves_with_all_pvals(
    gene_stats: pd.DataFrame,
    groupby_cols: list[str],
    perturbed_col: str = "perturbed",
    score_col: str = "-log10_pval",
) -> pd.DataFrame:
    curves = get_curves(gene_stats, groupby_cols, perturbed_col, score_col)
    return left_join_other_pvals_to_curves(curves, gene_stats, score_col)


def get_curves(
    gene_stats: pd.DataFrame,
    groupby_cols: list[str],
    perturbed_col: str = "perturbed",
    score_col: str = "-log10_pval",
) -> pd.DataFrame:
    def f(df: pd.DataFrame) -> pd.DataFrame:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fpr, tpr, scores = roc_curve(
                y_true=df[perturbed_col],
                y_score=df[score_col],
            )
        # assert they are the same shape
        assert fpr.shape == tpr.shape == scores.shape
        return pd.DataFrame(
            {
                "fpr": fpr,
                "tpr": tpr,
                "n": len(df),
                "n_true": df[perturbed_col].sum(),
                "n_false": len(df) - df[perturbed_col].sum(),
            },
            index=pd.Index(scores, name=score_col),
        )

    # compute curves #############
    dfg = gene_stats.groupby(groupby_cols)
    logger.debug("calculating curves of classifier metrics with %s", score_col)
    df = dfg.apply(f)
    df["fp"] = df["fpr"] * df["n_false"]
    df["tp"] = df["tpr"] * df["n_true"]
    df["precision"] = df["tp"] / (df["tp"] + df["fp"])
    df["recall"] = df["tpr"]
    # results in TypeError... can't convert NaN to <NA> for some reason
    # df = df.astype({"fp": "Int64", "tp": "Int64"})
    return df


def left_join_other_pvals_to_curves(
    curves: pd.DataFrame,
    gene_stats: pd.DataFrame,
    score_col: str = "-log10_pval",
):
    # get other pvals #############
    assert score_col in curves.index.names
    groupby_cols = list(curves.index.names)
    columns = gene_stats.columns
    other_pval_cols = list(filter(lambda x: "pval" in x and x != score_col, columns))
    logger.debug("other pvals: %s", other_pval_cols)
    df_other_pvals = gene_stats.groupby(groupby_cols)[other_pval_cols].first()

    # add other pvals with left join #############
    # assert index level names are the same
    assert curves.index.names == df_other_pvals.index.names, (
        curves.index.names,
        df_other_pvals.index.names,
    )
    # print shapes of each
    print(curves.shape, df_other_pvals.shape)
    logger.debug("adding other pvals: %s", other_pval_cols)
    return curves.join(df_other_pvals, how="left", validate="one_to_one")
