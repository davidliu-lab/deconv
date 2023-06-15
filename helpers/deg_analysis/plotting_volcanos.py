import logging

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from helpers.deg_analysis.plotting_utils import (
    add_fdr_lines,
    remove_excess_facet_axis_titles,
)

logger = logging.getLogger(__name__)


def make_volcano_grid_scatter(
    gene_stats: pd.DataFrame,
    aggregate_by: tuple[str],
    pval_col: str = "-log10_pval_adjusted_bh",
    perturbed_col: str = "gene_perturbed",
    marker_color: str | None = None,
    hover_data: tuple[str] = (),
    facet_col: str | None = "log2_fc",
    facet_row: str | None = "malignant_means",
    simplify_facet_titles: bool = True,
    with_fdr_lines: bool = True,
) -> go.Figure:
    assert facet_col in aggregate_by
    assert facet_row in aggregate_by
    dfg = gene_stats.groupby(aggregate_by)
    logger.debug("Distribution of group sizes %s", dfg.size().value_counts())
    marker_color_tuple = (marker_color,) if marker_color else ()
    aggregands = ("log2_fold_change", pval_col, perturbed_col) + hover_data + marker_color_tuple
    logger.debug("Aggregating %s by %s", aggregands, aggregate_by)
    df_to_plot = dfg[aggregands].median().reset_index()
    logger.debug("plotting df_to_plot with shape %s", df_to_plot.shape)
    if not marker_color:
        marker_color = perturbed_col
    fig = px.scatter(
        df_to_plot,
        x="log2_fold_change",
        y=pval_col,
        color=marker_color,
        symbol=perturbed_col,
        hover_name="gene_symbol",
        hover_data=hover_data,
        facet_row=facet_row,
        facet_col=facet_col,
        symbol_map={False: "circle", True: "x"},
    )
    fig = format_volcano_figure(fig, simplify_facet_titles, with_fdr_lines)
    return fig


def format_volcano_figure(
    fig: go.Figure,
    simplify_facet_titles: bool = True,
    with_fdr_lines: bool = True,
    marker_size: int = 5,
) -> go.Figure:
    try:
        fig.update_traces(marker_size=marker_size)
    except ValueError:
        pass
    fig.update_xaxes(range=[-2, 2])
    fig.update_yaxes(range=[0, 6])
    if simplify_facet_titles:
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
        fig = remove_excess_facet_axis_titles(fig)
    if with_fdr_lines:
        fig = add_fdr_lines(fig, axis="y")
    fig.update_layout(width=800, height=500)
    # overlay legend on top right
    fig.update_layout(
        legend=dict(yanchor="top", y=0.95, xanchor="right", x=0.95),
    )
    return fig
