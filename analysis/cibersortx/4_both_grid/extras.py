from functools import partial
import logging
import re

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import upath

import sklearn
import sklearn.metrics


logger = logging.getLogger(__name__)


def get_parquet_paths(path_root: upath.UPath) -> list[upath.UPath]:
    paths = path_root.glob(
        "**/deg_analysis/gene_stats_malignant_cibersortx.parquet"
    )
    paths = map(str, paths)
    # paths = filter(re.compile(r".*run_id=0[0-2].*").match, paths)
    paths = sorted(paths)
    paths = list(paths)
    return paths


def extract_from_path(path: str, var_name: str) -> str:
    _ = path.split(var_name + "=")[1]
    return _.split("/")[0]


def extract_vars_from_path(path: str) -> list[tuple[str, str]]:
    # extract the variables from the path
    return re.findall(r"(\w+)=(\w+)", path)


def test_extract_vars_from_path():
    test_path = "thing1=0/a=1,b=2/foo=bar/thing.parquet"
    result = extract_vars_from_path(test_path)
    expectation = [("thing1", "0"), ("a", "1"), ("b", "2"), ("foo", "bar")]
    assert result == expectation


orderings = {
    "malignant_means": [
        "0.71,0.71",
        "0.7,0.72",
        "0.65,0.75",
        "0.6,0.8",
        "0.55,0.85",
    ],
    "log2_fc": [
        "-1.50",
        "-1.00",
        "-0.50",
        "-0.25",
        "0.00",
        "0.25",
        "0.50",
        "1.00",
        "1.50",
    ],
}


def load_gene_stats(path_root: upath.UPath):
    parquet_paths = get_parquet_paths(path_root)
    df = pd.concat(
        {str(path): pd.read_parquet(path) for path in parquet_paths},
        names=["path", "index"],
    )

    for column_name in ["malignant_means", "log2_fc", "run_id"]:
        s = df.index.get_level_values("path").map(
            partial(extract_from_path, var_name=column_name)
        )
        if column_name in orderings:
            logger.debug("Setting ordering for %s", column_name)
            dtype = pd.CategoricalDtype(orderings[column_name], ordered=True)
            s = s.astype(dtype)
        df[column_name] = s

    # add gene_perturbed column
    path_genes_perturbed = (
        upath.UPath(
            "gs://liulab/differential_composition_and_expression/20230224_07h54m40s"
        )
        / "genes_perturbed.csv"
    )
    genes_perturbed = pd.read_csv(
        path_genes_perturbed,
        index_col=0,
    )
    df["gene_perturbed"] = df["gene_symbol"].isin(genes_perturbed.index)

    df = df.set_index(
        [
            "malignant_means",
            "log2_fc",
            "run_id",
            "gene_symbol",
            "gene_perturbed",
        ]
    )
    df = df.sort_index()
    return df


def make_volcano_grid_scatter(df: pd.DataFrame) -> go.Figure:
    fields_to_groupby = [
        "malignant_means",
        "log2_fc",
        "gene_symbol",
        "gene_perturbed",
    ]
    dfg = df.groupby(fields_to_groupby)
    fields_to_aggregate = [
        "log2_fold_change",
        "-log10_pval",
        "significant_bh_fdr=0.10",
    ]
    df = dfg[fields_to_aggregate].mean()
    df = df.reset_index()
    fig = px.scatter(
        df,
        x="log2_fold_change",
        y="-log10_pval",
        facet_col="malignant_means",
        facet_row="log2_fc",
        # hover_name="gene_symbol",
        color="gene_perturbed",
    )
    # set marker size to 1
    fig.update_traces(marker_size=1)
    fig.update_xaxes(range=[-2, 2])
    fig.update_yaxes(range=[0, 6])
    return fig


def make_volcano_grid_density_contour(
    df: pd.DataFrame,
) -> go.Figure:
    fields_to_groupby = [
        "malignant_means",
        "log2_fc",
        "gene_symbol",
        "gene_perturbed",
    ]
    dfg = df.groupby(fields_to_groupby)
    fields_to_aggregate = [
        "log2_fold_change",
        "-log10_pval",
        "significant_bh_fdr=0.10",
    ]
    df = dfg[fields_to_aggregate].mean()
    df = df.reset_index()
    fig = px.density_contour(
        df,
        x="log2_fold_change",
        y="-log10_pval",
        facet_col="malignant_means",
        facet_row="log2_fc",
        # hover_name="gene_symbol",
        # color="gene_perturbed",
    )
    fig.update_xaxes(range=[-2, 2])
    fig.update_yaxes(range=[0, 6])
    return fig


def calculate_precision_and_recall(df: pd.DataFrame) -> pd.DataFrame:
    """
    - Computes precision and recall using the column "-log10_pval" as the score.
    - Uses sklearn.metrics.precision_recall_curve
    """
    df = df.reset_index()

    def compute_curve(df: pd.DataFrame) -> pd.DataFrame:
        precision, recall, thresholds = sklearn.metrics.precision_recall_curve(
            y_true=df["gene_perturbed"],
            probas_pred=df["pval"],
        )
        # extend thresholds by one to include infinity
        thresholds = np.append(thresholds, np.inf)
        return pd.DataFrame(
            {
                "precision": precision,
                "recall": recall,
                "thresholds": thresholds,
            }
        )

    dfg = df.groupby(["malignant_means", "log2_fc", "run_id"])
    df = dfg.apply(compute_curve)
    return df


def calculate_roc(df: pd.DataFrame) -> pd.DataFrame:
    """
    - Computes precision and recall using the column "-log10_pval" as the score.
    - Uses sklearn.metrics.precision_recall_curve
    """
    df = df.reset_index()

    def compute_curve(df: pd.DataFrame) -> pd.DataFrame:
        fpr, tpr, thresholds = sklearn.metrics.roc_curve(
            y_true=df["gene_perturbed"],
            y_score=df["-log10_pval"],
        )
        # extend thresholds by one to include infinity
        # thresholds = np.append(thresholds, np.inf)
        logger.debug(
            "shapes: %s, %s, %s", fpr.shape, tpr.shape, thresholds.shape
        )
        return pd.DataFrame(
            {
                "fpr": fpr,
                "tpr": tpr,
                "thresholds": thresholds,
            }
        )

    dfg = df.groupby(["malignant_means", "log2_fc", "run_id"])
    df = dfg.apply(compute_curve)
    return df


def plot_roc(df: pd.DataFrame) -> go.Figure:
    df = df.reset_index()
    fig = px.line(
        df,
        x="fpr",
        y="tpr",
        labels={"x": "False Positive Rate", "y": "True Positive Rate"},
        facet_col="malignant_means",
        facet_row="log2_fc",
        hover_data=["thresholds"],
        color="run_id",
    )
    fig.add_shape(type="line", line=dict(dash="dash"), x0=0, x1=1, y0=0, y1=1)
    fig.update_xaxes(constrain="domain", dtick=0.2)
    fig.update_yaxes(scaleanchor="x", scaleratio=1, dtick=0.2)
    fig.update_xaxes(range=[0, 1])
    fig.update_yaxes(range=[0, 1])
    fig.update_layout(
        title="ROC Curve",
    )
    return fig


def plot_precision_recall_curve(df: pd.DataFrame) -> go.Figure:
    df = df.reset_index()
    fig = px.line(
        df,
        x="recall",
        y="precision",
        labels={"x": "Recall", "y": "Precision"},
        facet_col="malignant_means",
        facet_row="log2_fc",
        hover_data=["thresholds"],
        color="run_id",
    )
    fig.update_xaxes(range=[0, 1])
    fig.update_yaxes(range=[0, 1])
    fig.update_layout(
        title="Precision-Recall Curve",
    )
    return fig
