from functools import partial
import itertools
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
print("loading extras module")


def get_parquet_paths(path_root: upath.UPath) -> list[upath.UPath]:
    paths = path_root.glob("**/gene_stats_*")
    paths = map(str, paths)
    # paths = filter(re.compile(r".*run_id=0[0-2].*").match, paths)
    paths = sorted(paths)
    paths = list(paths)
    return paths


def extract_from_path(path: str, var_name: str) -> str:
    _ = path.split(var_name + "=")[1]
    return _.split("/")[0]


def test_extract_from_path():
    result = extract_from_path("thing/a=10", "a")
    assert result == "10"


def extract_vars_from_path(path: str) -> list[tuple[str, str]]:
    # extract the variables from the path
    return re.findall(r"(\w+)=(\w+)", path)


def test_extract_vars_from_path():
    test_path = "thing1=0/a=1,b=2/foo=bar/thing.parquet"
    result = extract_vars_from_path(test_path)
    expectation = [("thing1", "0"), ("a", "1"), ("b", "2"), ("foo", "bar")]
    assert result == expectation


def extract_origin_from_path(path: str) -> str:
    # extract the substring between "gene_stats_" and ".parquet"
    return re.findall(r"gene_stats_(.+)\.parquet", path)[0]


def test_extract_origin_from_path():
    test_path = "/foo=bar/gene_stats_thing.parquet"
    result = extract_origin_from_path(test_path)
    assert result == "thing"


ordering_functions = {
    "log2_fc": float,
    "malignant_means": lambda x: "0.715,0.715" if x == "None,None" else x,
}


def load_gene_stats(path_root: upath.UPath):
    parquet_paths = get_parquet_paths(path_root)
    logger.debug(parquet_paths)
    df = pd.concat(
        {str(path): pd.read_parquet(path) for path in parquet_paths},
        names=["path", "index"],
    )

    for column_name in ["malignant_means", "log2_fc", "run_id"]:
        try:
            s = df.index.get_level_values("path").map(
                partial(extract_from_path, var_name=column_name)
            )
        except IndexError:
            logger.debug("skipping %s because not found", column_name)
            continue
        if column_name in ordering_functions:
            logger.debug("Setting ordering for %s", column_name)
            values = s.unique()
            values = sorted(values, key=ordering_functions[column_name])
            dtype = pd.CategoricalDtype(values, ordered=True)
            s = s.astype(dtype)
        df[column_name] = s

    # add "origin" column
    df["origin"] = df.index.get_level_values("path").map(extract_origin_from_path)

    # add "gene_perturbed" column
    path_genes_perturbed = (
        upath.UPath("gs://liulab/differential_composition_and_expression/20230224_07h54m40s")
        / "genes_perturbed.csv"
    )
    genes_perturbed = pd.read_csv(
        path_genes_perturbed,
        index_col=0,
    )
    df["gene_perturbed"] = df["gene_symbol"].isin(genes_perturbed.index)
    # if log2_fc is 0, then gene_perturbed should be False
    df["gene_perturbed"] = df["gene_perturbed"] & df["log2_fc"].astype(float).abs() > 0

    df = df.set_index(
        [
            "origin",
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
        facet_col="log2_fc",
        facet_row="malignant_means",
        # hover_name="gene_symbol",
        color="gene_perturbed",
        # use marker symbols "." and "*" for "gene_perturbed" values of False and True
        # respectively
        symbol="gene_perturbed",
        symbol_map={False: "circle", True: "x"},
    )
    # set marker size to 1
    fig.update_traces(marker_size=1)
    fig.update_xaxes(range=[-2, 2])
    fig.update_yaxes(range=[0, 6])
    # remove variable name from facet label
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
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
        facet_col="log2_fc",
        facet_row="malignant_means",
        # hover_name="gene_symbol",
        # color="gene_perturbed",
    )
    fig.update_xaxes(range=[-2, 2])
    fig.update_yaxes(range=[0, 6])
    # remove variable name from facet label
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    return fig


def compute_classification_score(df: pd.DataFrame) -> pd.Series:
    sign_observed_log2_fc = np.sign(df["log2_fold_change"]).replace({0: 1})
    sign_dist_log2_fc = np.sign(df["log2_fc"].astype(float)).replace({0: 1})
    return df["-log10_pval"] * sign_dist_log2_fc * sign_observed_log2_fc


def calculate_roc(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    - Computes precision and recall using the column "-log10_pval" as the score.
    - Uses sklearn.metrics.roc_curve
    """
    df = df.reset_index()
    df["classification_score"] = compute_classification_score(df)

    def compute_curve_for_group(df: pd.DataFrame) -> pd.DataFrame:
        data = {}
        data["fpr"], data["tpr"], data["thresholds"] = sklearn.metrics.roc_curve(
            y_true=df["gene_perturbed"],
            y_score=df["classification_score"],
        )
        return pd.DataFrame.from_records(data)

    def compute_roc_auc_score_for_group(df: pd.DataFrame) -> float:
        try:
            return sklearn.metrics.roc_auc_score(
                y_true=df["gene_perturbed"],
                y_score=df["classification_score"],
            )
        except ValueError:
            return np.nan

    dfg = df.groupby(["malignant_means", "log2_fc", "run_id"])
    roc_curves = dfg.apply(compute_curve_for_group)
    roc_auc_scores = dfg.apply(compute_roc_auc_score_for_group)
    return roc_curves, roc_auc_scores


def plot_roc(df: pd.DataFrame) -> go.Figure:
    df = df.reset_index()
    fig = px.line(
        df,
        x="fpr",
        y="tpr",
        labels={"x": "False Positive Rate", "y": "True Positive Rate"},
        facet_col="log2_fc",
        facet_row="malignant_means",
        hover_data=["thresholds"],
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
    return fig


def calculate_precision_and_recall(df: pd.DataFrame) -> pd.DataFrame:
    """
    - Computes precision and recall using the column "-log10_pval" as the score.
    - Uses sklearn.metrics.precision_recall_curve
    """
    df = df.reset_index()
    df["classification_score"] = compute_classification_score(df)

    def compute_curve(df: pd.DataFrame) -> pd.DataFrame:
        precision, recall, thresholds = sklearn.metrics.precision_recall_curve(
            y_true=df["gene_perturbed"],
            probas_pred=df["classification_score"],
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


def plot_precision_recall_curve(df: pd.DataFrame) -> go.Figure:
    df = df.reset_index()
    fig = px.line(
        df,
        x="recall",
        y="precision",
        labels={"x": "Recall", "y": "Precision"},
        facet_col="log2_fc",
        facet_row="malignant_means",
        hover_data=["thresholds"],
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
        title="Precision-Recall Curve",
        showlegend=False,
    )
    # remove variable name from facet label
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    return fig
