import itertools
import logging

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from helpers.deg_analysis.plotting_utils import _util_remove_excess_axis_titles

logger = logging.getLogger(__name__)


def plot_roc(df: pd.DataFrame) -> go.Figure:
    fig = px.line(
        df.reset_index(),
        x="fpr",
        y="tpr",
        labels={"x": "False Positive Rate", "y": "True Positive Rate"},
        facet_col="log2_fc",
        facet_row="malignant_means",
        hover_data=list(df.columns),
        color="run_id",
    )
    fig.update_xaxes(
        range=[0, 1],
        constrain="domain",
        dtick=0.2,
    )
    fig.update_yaxes(
        range=[0, 1],
        scaleanchor="x",
        scaleratio=1,
        dtick=0.2,
    )
    fig.update_layout(
        title="ROC Curves",
        showlegend=False,
    )
    # remove variable name from facet label
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    # add dashed diagonal line in each plot
    for row, col in itertools.product(*fig._get_subplot_rows_columns()):
        fig.add_shape(type="line", line=dict(dash="dash"), x0=0, x1=1, y0=0, y1=1, row=row, col=col)
    fig.update_layout(width=800, height=800)
    return fig


def plot_precision_recall_curve(df: pd.DataFrame) -> go.Figure:
    # df = df.loc[
    #     df.index.get_level_values("malignant_means").isin(
    #         ["0.55,0.85", "0.65,0.75", "None,None", "0.75,0.65", "0.85,0.55"]
    #     )
    # ]
    # df = df.loc[df.index.get_level_values("log2_fc").isin(["-0.50", "0.50"])]
    # df = df.loc[df.index.get_level_values("run_id").isin(["00", "01"])]
    fig = px.line(
        df.reset_index(),
        x="recall",
        y="precision",
        # labels={"x": "Recall", "y": "Precision"},
        facet_col="log2_fc",
        facet_row="malignant_means",
        hover_data=list(df.columns),
        color="run_id",
    )
    fig.update_xaxes(
        range=[0, 1],
        constrain="domain",
        dtick=0.2,
    )
    fig.update_yaxes(
        range=[0, 1],
        scaleanchor="x",
        scaleratio=1,
        dtick=0.2,
    )
    fig.update_layout(
        title="precision & recall",
        showlegend=False,
    )
    # remove variable name from facet label
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    # fig = _util_remove_excess_axis_titles(fig)
    fig.update_layout(width=800, height=800)
    return fig


def plot_metric_by_threshold(
    df_curves: pd.DataFrame, score_column: str, metric_column: str
) -> go.Figure:
    fig = px.line(
        df_curves.reset_index(),
        x=score_column,
        y=metric_column,
        color="run_id",
        markers=True,
        facet_col="log2_fc",
        facet_row="malignant_means",
    )
    # add vertical line at x=1.0
    if score_column.startswith("-log10_pval_adjusted_bh"):
        for row, col in fig._get_subplot_coordinates():
            fig.add_shape(
                type="line",
                x0=1.0,
                y0=0,
                x1=1.0,
                y1=1,
                line=dict(color="black", width=2, dash="dot"),
                row=row,
                col=col,
            )
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    fig = _util_remove_excess_axis_titles(fig)
    return fig


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


def make_score_table_with_stddev(scores: pd.Series, cmap: str = "RdBu") -> pd.DataFrame:
    aggfuncs = ["mean", "std"]
    aggregations = scores.groupby(["malignant_means", "log2_fc"]).agg(func=aggfuncs)
    means_and_stddevs = aggregations.apply(
        lambda row: f"{row['mean']:4.2f}±{row['std']:4.2f}", axis="columns"
    )
    df_means = aggregations["mean"].unstack("log2_fc")
    df_means_and_stddevs = means_and_stddevs.unstack("log2_fc")
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
        df_scores.groupby(["malignant_means", "log2_fc"])
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
