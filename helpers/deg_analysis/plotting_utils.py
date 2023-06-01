import logging
import statistics

import plotly.graph_objects as go

logger = logging.getLogger(__name__)


def _util_remove_excess_axis_titles(fig: go.Figure) -> go.Figure:
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
