import logging
import statistics

import numpy as np
import plotly.graph_objects as go

logger = logging.getLogger(__name__)


def remove_excess_facet_axis_titles(fig: go.Figure) -> go.Figure:
    # remove axis titles from all subplots... except the middle one

    coordinates = list(fig._get_subplot_coordinates())
    middle_row = statistics.median_low([row for row, _ in coordinates])
    middle_col = statistics.median_low([col for _, col in coordinates])
    logger.debug("middle_col=%s", middle_col)
    logger.debug("middle_row=%s", middle_row)
    for row, col in coordinates:
        if row != middle_row:
            fig.update_yaxes(title_text="", row=row, col=col)
        if col != middle_col:
            fig.update_xaxes(title_text="", row=row, col=col)
    return fig


def remove_variable_names_from_facet_axis_titles(fig: go.Figure) -> go.Figure:
    return fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))


def add_fdr_lines(
    fig: go.Figure,
    alpha: float = 0.1,
    log_scale: bool = True,
    axis: str = "x",
) -> go.Figure:
    """Add FDR lines to a figure"""
    if log_scale:
        value = -np.log10(alpha)
    else:
        value = alpha
    for row, col in fig._get_subplot_coordinates():
        # annotation is float with no trailing zeros
        kwargs = dict(
            line_dash="dot",
            annotation_text=f"FDR={alpha}",
            row=row,
            col=col,
        )
        if axis == "x":
            fig.add_vline(x=value, **kwargs)
        elif axis == "y":
            fig.add_hline(y=value, **kwargs)
    return fig
