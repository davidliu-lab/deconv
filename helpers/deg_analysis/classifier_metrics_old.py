"""
Old functions for computing classifier metrics (ROC, PR, etc.) from a dataframe of
gene-level statistics.
"""
import logging
import warnings

import numpy as np
import pandas as pd
import sklearn
import sklearn.metrics
from pandas.core.groupby.generic import DataFrameGroupBy
from sklearn.metrics import (
    precision_recall_curve,
    precision_recall_fscore_support,
    precision_score,
    roc_auc_score,
    roc_curve,
)

logger = logging.getLogger(__name__)


def _get_groupby(df: pd.DataFrame) -> DataFrameGroupBy:
    # note - this is extremely fast. no need to reuse the result.
    if "gene_symbol" in df.index.names:
        groupby_fields = list(filter(lambda x: x != "gene_symbol", df.index.names))
    else:
        groupby_fields = ["origin", "malignant_means", "log2_fc", "run_id"]
    logger.debug("grouping by %s", groupby_fields)
    return df.groupby(
        groupby_fields,
        # these just remove all index levels, don't include index defined inside groupby.apply
        # group_keys=True,
        # as_index=False,
    )


def calculate_all_curves(
    df: pd.DataFrame,
    score_col: str,
    perturbed_col: str = "gene_perturbed",
) -> pd.DataFrame:
    """Calculate fpr, tpr, fnr, tnr curves for each group in df"""
    assert df[perturbed_col].dtype == bool, df[perturbed_col].dtype
    dfg = _get_groupby(df)

    def compute_curves_for_group(df: pd.DataFrame) -> pd.DataFrame:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fpr, tpr, scores = sklearn.metrics.roc_curve(
                y_true=df[perturbed_col],
                y_score=df[score_col],
            )
            n_actual_positives = df[perturbed_col].sum()
            if n_actual_positives == 0:
                tp = np.zeros_like(tpr)
                fn = np.zeros_like(tpr)
            else:
                tp = tpr * n_actual_positives
                fn = n_actual_positives - tp
            fp = fpr * (len(df) - n_actual_positives)
            tn = len(df) - n_actual_positives - fp
            if n_actual_positives != 0:
                assert np.allclose(tp + fn, n_actual_positives), tp + fn
                assert np.allclose(tpr, tp / (tp + fn))
            assert np.allclose(fp + tn, len(df) - n_actual_positives), fp + tn
            # assert no NaNs, infinities, etc.
            assert np.all(np.isfinite(tp))
            assert np.all(np.isfinite(fp))
            assert np.all(np.isfinite(fn))
            assert np.all(np.isfinite(tn))
            result = pd.DataFrame(
                {
                    # score_col: scores,
                    "fpr": fpr,
                    "tpr": tpr,
                    "precision": tp / (tp + fp),
                    "recall": tpr,
                    "fp": fp,
                    "tp": tp,
                    "fn": fn,
                    "tn": tn,
                },
                index=pd.Index(scores, name=score_col),
            )
        # result = result.rename_axis(index="ordering")
        # result = result.set_index(score_col)
        return result

    logger.debug("calculating curves of classifier metrics with %s", score_col)
    df_result = dfg.apply(compute_curves_for_group)
    df_result = df_result.astype({"fp": int, "tp": int, "fn": int, "tn": int})
    return df_result


def calculate_roc_curves(
    df: pd.DataFrame,
    score_col: str,
    perturbed_col: str = "gene_perturbed",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    dfg = _get_groupby(df)

    def compute_roc_curve_for_group(df: pd.DataFrame) -> pd.DataFrame:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fpr, tpr, scores = roc_curve(
                y_true=df[perturbed_col],
                y_score=df[score_col],
            )
        return pd.DataFrame({"fpr": fpr, "tpr": tpr, score_col: scores})

    def compute_roc_auc_score_for_group(df: pd.DataFrame) -> float:
        try:
            return roc_auc_score(
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


def calculate_precision_and_recall_curves(
    df: pd.DataFrame,
    score_col: str,
    perturbed_col: str = "gene_perturbed",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    dfg = _get_groupby(df)

    def compute_precision_recall_curve_for_group(df: pd.DataFrame) -> pd.DataFrame:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            precision, recall, scores = precision_recall_curve(
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

    def compute_precision_for_group(df: pd.DataFrame) -> float:
        # with warnings.catch_warnings():
        #     warnings.simplefilter("ignore")
        return precision_score(
            y_true=df[perturbed_col],
            y_pred=df["significant_bh_fdr=0.10"],
            zero_division=0,
        )

    logger.debug("calculating precision-recall curves with %s", score_col)
    df_curves = dfg.apply(compute_precision_recall_curve_for_group)
    logger.debug("calculating precision values with %s", score_col)
    df_precision = dfg.apply(compute_precision_for_group)
    return df_curves, df_precision


def compute_scores(
    df: pd.DataFrame | pd.core.groupby.generic.DataFrameGroupBy,
    perturbed_col: str,
) -> pd.DataFrame:
    """Compute precision, recall and ROC AUC scores for each experimental group"""

    def f(df: pd.DataFrame) -> pd.Series:
        y_true = df[perturbed_col]
        y_pred = df["significant_bh_fdr=0.10"]
        y_score = df["-log10_pval"]
        precision, recall, _, _ = precision_recall_fscore_support(
            y_true=y_true,
            y_pred=y_pred,
            zero_division=0,
            average="binary",
        )
        try:
            roc_auc = roc_auc_score(
                y_true=y_true,
                y_score=y_score,
            )
        except ValueError:
            roc_auc = np.nan
        data = {
            "roc_auc": roc_auc,
            "precision": precision,
            "recall": recall,
            "tp": np.sum(y_true & y_pred),
            "fp": np.sum(~y_true & y_pred),
            "tn": np.sum(~y_true & ~y_pred),
            "fn": np.sum(y_true & ~y_pred),
        }
        return pd.Series(data)

    if isinstance(df, pd.DataFrame):
        dfg = _get_groupby(df)
    else:
        dfg = df
    logger.debug("computing scores")
    df_scores = dfg.apply(f)
    logger.debug("converting counts to int")
    count_cols = list(filter(lambda x: x.endswith("_count"), df_scores.columns))
    df_scores[count_cols] = df_scores[count_cols].astype(int)
    return df_scores


def compute_scores_for_threshold(
    gene_stats: pd.DataFrame | pd.core.groupby.generic.DataFrameGroupBy,
    perturbed_col: str,
    score_col: str,
    threshold: float,
) -> pd.DataFrame:
    """Compute precision, recall and ROC AUC scores for each experimental group"""

    def f(df: pd.DataFrame) -> pd.Series:
        y_true = df[perturbed_col]
        y_pred = df[score_col] >= threshold
        data = {}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            precision, recall, _, _ = precision_recall_fscore_support(
                y_true=y_true,
                y_pred=y_pred,
                # zero_division=0,  # why not np.nan? can only be 0, 1, or "warn"
                average="binary",
            )
        data["precision"], data["recall"] = precision, recall
        data["tp"] = tp = np.sum(y_true & y_pred)
        data["fp"] = fp = np.sum(~y_true & y_pred)
        data["tn"] = tn = np.sum(~y_true & ~y_pred)
        data["fn"] = fn = np.sum(y_true & ~y_pred)
        if np.sum(y_true) == 0:
            assert precision == 0.0
            assert recall == 0.0
            data["fpr"] = 0.0
        else:
            assert precision == tp / (
                tp + fp
            ), f"{df.name}, {precision} != {tp} / ({tp} + {fp}), {np.sum(y_true)}"
            assert recall == tp / (tp + fn), f"{recall} != {tp} / ({tp} + {fn})"
            data["fpr"] = fp / (fp + tn)
        data["tpr"] = recall
        data["tnr"] = 1.0 - data["fpr"]
        data["fnr"] = 1.0 - data["tpr"]
        s = pd.Series(data)
        print(s)
        return s

    if isinstance(gene_stats, pd.DataFrame):
        dfg = _get_groupby(gene_stats)
    else:
        dfg = gene_stats
    logger.debug("computing scores")
    df_scores = dfg.apply(f)
    # logger.debug("converting counts to int")
    # df_scores = df_scores.astype({"tp": int, "fp": int, "tn": int, "fn": int})
    return df_scores


def compute_all_curves_and_metrics(
    df: pd.DataFrame, signed_directional: bool = True
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if signed_directional:
        score_column_roc = "-log10_pval_signed_directional"
        score_column_precision = "-log10_pval_adjusted_bh_signed_directional"
    else:
        score_column_roc = "-log10_pval"
        score_column_precision = "-log10_pval_adjusted_bh"

    df_roc_curves, _ = calculate_roc_curves(
        df,
        score_column_roc,
        perturbed_col="perturbed",
    )
    df_precision_recall_curves, _ = calculate_precision_and_recall_curves(
        df,
        score_column_precision,
        perturbed_col="perturbed",
    )
    df_scores = compute_scores(df, perturbed_col="perturbed")

    return (
        df_roc_curves,
        df_precision_recall_curves,
        df_scores,
    )


def get_metrics_at_threshold(
    df_curves: pd.DataFrame,
    groupby: list[str],
    score: str,
    threshold: float,
    metrics: list[str],
) -> pd.DataFrame:
    try:
        df_curves = df_curves.reset_index(level=score)
    except KeyError:
        pass
    # df_curves = df_curves[df_curves[score] >= threshold]  # this excludes runs with no significant genes!
    dfg = df_curves.groupby(groupby)

    def f(df: pd.DataFrame) -> pd.Series:
        # df = df.sort_values(score, ascending=False)  # significant genes first
        # not_significant = df[df[score] < threshold]  # exclude significant genes
        # # return the most significant entry under the threshold
        # return not_significant.iloc[0][metrics]
        df = df.sort_values(score)
        try:
            point = df[df[score] >= threshold].iloc[0]
        except IndexError:
            point = df.iloc[-1]
        return point[metrics]

    result = dfg.apply(f)
    return result


# bad example
# produces results that don't match manual calculation
"""
logger.info("computing metrics")
alphas = [
    # 0.05,
    0.1,
    # 0.25,
]

results = {
    alpha: classifier_metrics_old.get_metrics_at_threshold(
        df_curves=df_curves,
        groupby=groupby_cols,  # ["origin", "malignant_means", "log2_fc", "run_id"],
        score="-log10_pval_adjusted_bh_signed_directional",
        threshold=-1.0 * np.log10(alpha),
        # score="pval_adjusted_bh_signed_directional",  # doesn't work, because ordered in reverse
        # threshold=alpha,
        metrics=[
            "n",
            "n_true",
            "n_false",
            "fp",
            "fpr",
            "tp",
            "tpr",
            "precision",
            "recall",
            "pval_adjusted_bh_signed_directional",
            "-log10_pval_adjusted_bh_signed_directional",
        ],
    )
    for alpha in alphas
}
results = pd.concat(results, names=["alpha"])
# reorder index levels
results = results.reorder_levels(["malignant_means", "run_id", "alpha"])
results = results.sort_index()
results["-log10(alpha)"] = -1.0 * np.log10(results.index.get_level_values("alpha"))
# add as index level
results = results.set_index("-log10(alpha)", append=True)
results
"""
