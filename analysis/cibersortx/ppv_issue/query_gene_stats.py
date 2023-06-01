import logging

import duckdb
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from results import gene_stats

logger = logging.getLogger(__name__)
logger.debug("loading module %s", __name__)

sql = """
select
    origin,
    malignant_means,
    log2_fc,
    run_id,
    gene_symbol,
    perturbed,
    "pval",
    "-log10_pval",
    "pval_adjusted_bh",
    -1.0 * log10("pval_adjusted_bh") as "-log10_pval_adjusted_bh",
from gene_stats
where
    origin = 'malignant_cibersortx'
    and malignant_means = 'None,None'
    and log2_fc in ('-1.50', '-1.00')
    --and run_id = '0'
    and gene_symbol like 'HLA%'
--using sample 1%
;
"""


# precision and recall curves
from sklearn.metrics import precision_recall_curve


def compute_all_precision_recall_curves(
    df: pd.DataFrame, groupbys: list[str], threshold_field: str
) -> pd.DataFrame:
    def agg_func(df_group: pd.DataFrame) -> pd.DataFrame:
        logger.info("computing precision-recall curve for group=%s", df_group.name)
        p, r, t = precision_recall_curve(df_group["perturbed"], df_group[threshold_field])
        logger.debug(
            "p[0]=%s, p[-1]=%s, r[0]=%s, r[-1]=%s, t[0]=%s, t[-1]=%s",
            p[0],
            p[-1],
            r[0],
            r[-1],
            t[0],
            t[-1],
        )
        return pd.DataFrame(
            {
                "precision": p[:-1],
                "recall": r[:-1],
                threshold_field: t,
            }
        )

    # agg_func = lambda df: compute_precision_recall_curve(df, threshold_field)
    df_precision_recall = df.groupby(groupbys).apply(agg_func).reset_index()
    logger.debug("df_precision_recall.shape: %s", df_precision_recall.shape)
    return df_precision_recall


# px.line with precision and recall on y axis, threshold on x axis


def plot_precision_recall_by_threshold(
    df_precision_recall: pd.DataFrame, threshold_field: str
) -> go.Figure:
    assert threshold_field in df_precision_recall.columns
    assert "run_id" in df_precision_recall.columns
    assert "log2_fc" in df_precision_recall.columns
    fig = px.line(
        df_precision_recall,
        x=threshold_field,
        y=["precision", "recall"],
        color="run_id",
        facet_col="log2_fc",
        title=f"Precision and Recall vs. {threshold_field}",
        labels={
            "value": "Precision or Recall",
            "variable": "Metric",
            threshold_field: threshold_field,
        },
    )

    if threshold_field == "-log10_pval_adjusted_bh":
        logger.debug("adding vertical line at -log10(0.1) for FDR=0.1")
        fig.add_vline(x=1, line_dash="dash", annotation_text="FDR=0.1")
        # x limit (0, 5)
        fig.update_xaxes(range=[0, 5])

    fig.update_yaxes(range=[0, 1])
    fig.update_layout(width=1000, height=600)
    # legend on top
    fig.update_layout(legend=dict(yanchor="top", y=0.99, xanchor="right", x=0.99))
    return fig


def plot_precision_recall_curves(
    df_precision_recall: pd.DataFrame, threshold_field: str
) -> go.Figure:
    fig = px.line(
        df_precision_recall,
        x="recall",
        y="precision",
        color="run_id",
        hover_data=[threshold_field],
        title="Precision-Recall Curve",
        facet_col="log2_fc",
    )
    # 800 width, 600 height
    fig.update_layout(width=800, height=600)
    fig.update_yaxes(range=[0, 1])
    return fig


if __name__ == "__main__":
    logger.setLevel("DEBUG")
    logger.debug("running main")

    df = duckdb.sql(sql).df()
    logger.debug("df.shape: %s", df.shape)
    print(df.head())

    threshold_field = "-log10_pval_adjusted_bh"

    df_precision_recall = compute_all_precision_recall_curves(df, threshold_field)

    fig = plot_precision_recall_curves(df_precision_recall, threshold_field)
    # write to PNG with scale to 2x
    fig.write_image("precision_recall_curve.png", scale=2)

    fig = plot_precision_recall_by_threshold(df_precision_recall, threshold_field)
    # write to PNG with scale to 2x
    fig.write_image("precision_recall_by_threshold.png", scale=2)
