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


def make_table_scores(df_scores: pd.DataFrame) -> pd.DataFrame:
    return (
        df_scores.groupby(["malignant_means", "log2_fc"]).mean().unstack("log2_fc")
    ).style.background_gradient(cmap="RdBu", axis=None, vmin=-0.5, vmax=1.5)


def make_table_scores_with_stddev(df_scores: pd.DataFrame) -> pd.DataFrame:
    def mean_and_stddev(series: pd.Series):
        mean = series.mean()
        stddev = series.std()
        return f"{mean:4.2f}±{stddev:4.2f}"

    return (
        df_scores.groupby(["malignant_means", "log2_fc"])
        .agg(func=mean_and_stddev)
        .unstack("log2_fc")
    ).style.background_gradient(cmap="Blues", axis=None, vmin=0.0, vmax=1.0)
