import logging

import numpy as np
import pandas as pd
import sklearn
import sklearn.metrics

logger = logging.getLogger(__name__)


def calculate_roc(df: pd.DataFrame, score_column: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    # df = df.loc[df.index.get_level_values("log2_fc") != "0.00"]
    # df = df.loc[df.index.get_level_values("run_id").isin(["00", "01"])]
    every_index_level_except_gene_symbol = list(
        filter(lambda x: x != "gene_symbol", df.index.names)
    )
    logger.debug("grouping by %s", every_index_level_except_gene_symbol)
    dfg = df.groupby(every_index_level_except_gene_symbol)

    def compute_roc_curve_for_group(df: pd.DataFrame) -> pd.DataFrame:
        fpr, tpr, scores = sklearn.metrics.roc_curve(
            y_true=df["gene_perturbed"],
            y_score=df[score_column],
        )
        return pd.DataFrame({"fpr": fpr, "tpr": tpr, score_column: scores})

    def compute_roc_auc_score_for_group(df: pd.DataFrame) -> float:
        try:
            return sklearn.metrics.roc_auc_score(
                y_true=df["gene_perturbed"],
                y_score=df[score_column],
            )
        except ValueError:
            return np.nan

    roc_curves = dfg.apply(compute_roc_curve_for_group)
    roc_auc_scores = dfg.apply(compute_roc_auc_score_for_group)
    return roc_curves, roc_auc_scores


def calculate_precision_and_recall(df: pd.DataFrame, score_column: str) -> pd.DataFrame:
    # df = df.loc[df.index.get_level_values("log2_fc") != "0.00"]
    every_index_level_except_gene_symbol = list(
        filter(lambda x: x != "gene_symbol", df.index.names)
    )
    logger.debug("grouping by %s", every_index_level_except_gene_symbol)
    dfg = df.groupby(every_index_level_except_gene_symbol)

    def compute_precision_recall_curve_for_group(df: pd.DataFrame) -> pd.DataFrame:
        precision, recall, scores = sklearn.metrics.precision_recall_curve(
            y_true=df["gene_perturbed"],
            probas_pred=df[score_column],
        )
        # extend scores by one to include infinity
        scores = np.append(scores, np.inf)
        return pd.DataFrame(
            {
                "precision": precision,
                "recall": recall,
                score_column: scores,
            }
        )

    def compute_precision_for_group(df: pd.DataFrame) -> pd.DataFrame:
        return sklearn.metrics.precision_score(
            y_true=df["gene_perturbed"],
            y_pred=df["significant_bh_fdr=0.10"],
        )

    df_curves = dfg.apply(compute_precision_recall_curve_for_group)
    df_precision = dfg.apply(compute_precision_for_group)
    return df_curves, df_precision
