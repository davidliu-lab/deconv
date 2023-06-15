import logging

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from numpy import add

from helpers.deg_analysis import plotting_utils
from helpers.deg_analysis.plotting_utils import remove_excess_facet_axis_titles

logger = logging.getLogger(__name__)


def plot_curves(
    df_curves: pd.DataFrame,
    x: str = "fpr",
    y: str = "tpr",
    facet_col: str = "log2_fc",
    facet_row: str = "malignant_means",
    color: str = "run_id",
    title: str = "ROC curves",
    simplify_facet_titles: bool = True,
):
    fig = px.line(
        df_curves,
        x=x,
        y=y,
        # labels={"x": "False Positive Rate", "y": "True Positive Rate"},
        facet_col=facet_col,
        facet_row=facet_row,
        hover_data=list(df_curves.columns),
        color=color,
    )
    fig = format_curves_fig(fig)
    fig.update_layout(title=title)
    if simplify_facet_titles:
        fig = plotting_utils.remove_variable_names_from_facet_axis_titles(fig)
    return fig


def plot_precision_recall_curves(df_curves: pd.DataFrame) -> go.Figure:
    return plot_curves(
        df_curves,
        x="recall",
        y="precision",
        facet_col="log2_fc",
        facet_row="malignant_means",
        color="run_id",
        title="Precision-recall curves",
    )


def plot_roc_curves(df_curves: pd.DataFrame) -> go.Figure:
    fig = plot_curves(df_curves)
    fig = add_roc_dashed_lines(fig)
    return fig


plot_roc = plot_roc_curves


def add_roc_dashed_lines(fig: go.Figure) -> go.Figure:
    # add dashed diagonal line in each plot
    # for row, col in itertools.product(*fig._get_subplot_rows_columns()):  # following should work?
    for row, col in fig._get_subplot_coordinates():
        fig.add_shape(type="line", line=dict(dash="dash"), x0=0, x1=1, y0=0, y1=1, row=row, col=col)
    return fig


def format_curves_fig(fig: go.Figure) -> go.Figure:
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
        showlegend=False,
    )
    width, height = [200 * max(x) for x in fig._get_subplot_rows_columns()]
    fig.update_layout(width=width, height=height)
    return fig


def plot_metric_by_threshold(
    df_curves: pd.DataFrame,
    score_column: str,
    metric_column: str,
    simplify_facet_titles: bool = True,
    facet_col: str = "log2_fc",
    facet_row: str = "malignant_means",
) -> go.Figure:
    fig = px.line(
        df_curves.reset_index(),
        x=score_column,
        y=metric_column,
        color="run_id",
        markers=True,
        facet_col=facet_col,
        facet_row=facet_row,
    )
    fig = format_metric_by_threshold_figure(
        fig,
        add_fdr_lines=score_column.startswith("-log10_pval_adjusted_bh"),
        simplify_facet_titles=simplify_facet_titles,
    )
    return fig


def format_metric_by_threshold_figure(
    fig: go.Figure,
    add_fdr_lines: bool = True,
    simplify_facet_titles: bool = True,
) -> go.Figure:
    # add vertical line at FDR alpha=0.1
    if add_fdr_lines:
        fig.update_xaxes(range=[-1, 6])
        plotting_utils.add_fdr_lines(fig, alpha=0.1, axis="x")
    if simplify_facet_titles:
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
        fig = remove_excess_facet_axis_titles(fig)
    fig.update_layout(
        legend=dict(yanchor="top", y=0.95, xanchor="right", x=0.95),
    )
    fig.update_yaxes(range=[0, 1])
    fig.update_layout(width=800, height=400)
    return fig
