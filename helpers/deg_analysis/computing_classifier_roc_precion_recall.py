import logging
import warnings

import numpy as np
import pandas as pd
from pandas.core.groupby.generic import DataFrameGroupBy
import sklearn
import sklearn.metrics

logger = logging.getLogger(__name__)


def _get_groupby(df: pd.DataFrame) -> DataFrameGroupBy:
    if "gene_symbol" in df.index.names:
        groupby_fields = list(filter(lambda x: x != "gene_symbol", df.index.names))
    else:
        groupby_fields = ["origin", "malignant_means", "log2_fc", "run_id"]
    logger.debug("grouping by %s", groupby_fields)
    return df.groupby(groupby_fields)


def calculate_roc(
    df: pd.DataFrame,
    score_col: str,
    perturbed_col: str = "gene_perturbed",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    dfg = _get_groupby(df)

    def compute_roc_curve_for_group(df: pd.DataFrame) -> pd.DataFrame:
        fpr, tpr, scores = sklearn.metrics.roc_curve(
            y_true=df[perturbed_col],
            y_score=df[score_col],
        )
        return pd.DataFrame({"fpr": fpr, "tpr": tpr, score_col: scores})

    def compute_roc_auc_score_for_group(df: pd.DataFrame) -> float:
        try:
            return sklearn.metrics.roc_auc_score(
                y_true=df[perturbed_col],
                y_score=df[score_col],
            )
        except ValueError:
            return np.nan

    logger.debug("calculating ROC curves with %s", score_col)
    roc_curves = dfg.apply(compute_roc_curve_for_group)
    logger.debug("calculating ROC AUC scores with %s", score_col)
    roc_auc_scores = dfg.apply(compute_roc_auc_score_for_group)
    return roc_curves, roc_auc_scores


def calculate_precision_and_recall(
    df: pd.DataFrame,
    score_col: str,
    perturbed_col: str = "gene_perturbed",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    dfg = _get_groupby(df)

    def compute_precision_recall_curve_for_group(df: pd.DataFrame) -> pd.DataFrame:
        precision, recall, scores = sklearn.metrics.precision_recall_curve(
            y_true=df[perturbed_col],
            probas_pred=df[score_col],
        )
        # extend scores by one to include infinity
        scores = np.append(scores, np.inf)
        return pd.DataFrame(
            {
                "precision": precision,
                "recall": recall,
                score_col: scores,
            }
        )

    def compute_precision_for_group(df: pd.DataFrame) -> pd.DataFrame:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return sklearn.metrics.precision_score(
                y_true=df[perturbed_col],
                y_pred=df["significant_bh_fdr=0.10"],
            )

    logger.debug("calculating precision-recall curves with %s", score_col)
    df_curves = dfg.apply(compute_precision_recall_curve_for_group)
    logger.debug("calculating precision values with %s", score_col)
    df_precision = dfg.apply(compute_precision_for_group)
    return df_curves, df_precision


def compute_all_curves_and_metrics(
    df: pd.DataFrame, signed_directional: bool = True
) -> tuple[pd.DataFrame, ...]:
    if signed_directional:
        score_column_roc = "-log10_pval_signed_directional"
        score_column_precision = "-log10_pval_adjusted_bh_signed_directional"
    else:
        score_column_roc = "-log10_pval"
        score_column_precision = "-log10_pval_adjusted_bh"

    df_roc_curves, df_roc_auc_scores = calculate_roc(
        df,
        score_column_roc,
        perturbed_col="perturbed",
    )
    df_precision_recall_curves, df_precision = calculate_precision_and_recall(
        df,
        score_column_precision,
        perturbed_col="perturbed",
    )

    return (
        df_roc_curves,
        df_roc_auc_scores,
        df_precision_recall_curves,
        df_precision,
    )
