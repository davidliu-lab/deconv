import logging

import pandas as pd

logger = logging.getLogger(__name__)


def make_score_table(scores: pd.Series, cmap: str = "RdBu") -> pd.DataFrame:
    # if some float type
    if scores.dtype == float:
        cmap = cmap or "RdBu"
        return (
            scores.groupby(["malignant_means", "log2_fc"]).mean().unstack("log2_fc")
        ).style.background_gradient(cmap=cmap, axis=None, vmin=-0.5, vmax=1.5)
    # if int
    elif scores.dtype == int:
        cmap = cmap or "Blues"
        df = (
            scores.groupby(["malignant_means", "log2_fc"]).mean().unstack("log2_fc")
        ).style.background_gradient(cmap=cmap, axis=None, vmin=0, vmax=500)
        df = df.style.format("{:.1f}")
        return df


def make_score_table_with_stddev(
    values: pd.Series,
    cmap: str = "RdBu",
    index: str = "malignant_means",
    columns: str = "log2_fc",
) -> pd.DataFrame:
    aggfuncs = ["mean", "std"]
    aggregations = values.groupby([index, columns]).agg(func=aggfuncs)
    means_and_stddevs = aggregations.apply(
        lambda row: f"{row['mean']:4.2f}±{row['std']:4.2f}", axis="columns"
    )
    df_means = aggregations["mean"].unstack(columns)
    df_means_and_stddevs = means_and_stddevs.unstack(columns)
    return df_means_and_stddevs.style.background_gradient(
        cmap=cmap, axis=None, vmin=-0.5, vmax=1.5, gmap=df_means
    )


def make_score_table_with_stddev_1(scores: pd.Series, cmap: str = "RdBu") -> pd.DataFrame:
    def mean_and_stddev(series: pd.Series):
        return f"{series.mean():4.2f}±{series.std():4.2f}"

    aggfuncs = {
        "mean": "mean",
        "mean_and_stddev": mean_and_stddev,
    }
    aggregations = (
        scores.groupby(["malignant_means", "log2_fc"])
        .agg(func=aggfuncs)
        .rename_axis(columns="aggregation")
        .unstack("log2_fc")
    )
    means = aggregations.xs("mean", level=1, axis="columns")
    means_and_stddevs = aggregations.xs("mean_and_stddev", level=1, axis="columns")
    return means_and_stddevs.style.background_gradient(
        cmap=cmap,
        axis=None,
        vmin=0.0,
        vmax=1.0,
        gmap=means,
    )
