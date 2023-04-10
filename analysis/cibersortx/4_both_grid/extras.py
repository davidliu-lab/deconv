from cgitb import reset
import itertools
import logging
import re
from functools import partial
from turtle import width

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import sklearn
import sklearn.metrics
import upath
from statsmodels.stats.multitest import multipletests

logger = logging.getLogger(__name__)
print("loading extras module")


def get_parquet_paths(path_root: upath.UPath) -> list[upath.UPath]:
    paths = path_root.glob("**/gene_stats_*")
    paths = map(str, paths)
    # paths = filter(re.compile(r".*run_id=0[0-1].*").match, paths)
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
    df = pd.concat(
        {str(path): pd.read_parquet(path) for path in parquet_paths},
        names=["path", "index"],
    )

    # add columns from path (e.g. "/foo=bar")
    for column_name in ["malignant_means", "log2_fc", "run_id"]:
        try:
            s = df.index.get_level_values("path").map(
                partial(extract_from_path, var_name=column_name)
            )
        except IndexError:
            logger.warning("skipping %s because not found", column_name)
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
    genes_perturbed = pd.read_csv(
        path_root / "genes_perturbed.csv",
        index_col=0,
    )
    gene_in_perturbation_list = df["gene_symbol"].isin(genes_perturbed.index)
    gene_perturbation_was_nonzero = df["log2_fc"].astype(float).abs() > 0
    df["gene_perturbed"] = gene_in_perturbation_list & gene_perturbation_was_nonzero

    # add signed p-value fields
    sign_observed_log2_fc = np.sign(df["log2_fold_change"].fillna(0)).replace({0: 1})
    sign_dist_log2_fc = np.sign(df["log2_fc"].astype(float)).replace({0: 1})
    # df["-log10_pval_signed"] already exists
    # df["-log10_pval_signed"] = df["-log10_pval"] * sign_observed_log2_fc
    df["-log10_pval_signed_directional"] = (
        df["-log10_pval"] * sign_dist_log2_fc * sign_observed_log2_fc
    )

    df = df.set_index(
        [
            "origin",
            "malignant_means",
            "log2_fc",
            "run_id",
            "gene_symbol",
            # "gene_perturbed",  # why should this be in the index?
        ]
    )
    df = df.sort_index()

    # add BH-adjusted p-value (aka q-value) and derived fields
    df_pval_adjusted = compute_pval_adjusted_fields(df)
    df = add_pval_adjusted_fields(df, df_pval_adjusted)
    return df


def check_gene_stats(df_gene_stats: pd.DataFrame):
    # groupby_fields = "origin	malignant_means	log2_fc	run_id".split("\t")
    groupby_fields = list(filter(lambda x: x != "gene_symbol", df_gene_stats.index.names))
    df_gene_stats_groupby = df_gene_stats.groupby(groupby_fields)
    print(df_gene_stats_groupby.size())
    # grouping by origin, count number of unique values of gene_symbol
    print(
        df_gene_stats.groupby("origin").apply(
            lambda df: df.index.get_level_values("gene_symbol").nunique()
        )
    )
    for level in df_gene_stats.index.names:
        # print level and nunique
        print(level, df_gene_stats.index.get_level_values(level).nunique())
    return df_gene_stats


def make_volcano_grid_scatter(df: pd.DataFrame) -> go.Figure:
    fields_to_groupby = list(filter(lambda x: x != "run_id", df.index.names))
    fields_to_groupby.append("gene_perturbed")
    logger.debug("Grouping by %s", fields_to_groupby)
    dfg = df.groupby(fields_to_groupby)
    logger.debug("Groups have sizes %s", dfg.size())
    fields_to_aggregate = [
        "log2_fold_change",
        "-log10_pval",
        "-log10_pval_adjusted_bh",
        "significant_bh_fdr=0.10",
    ]
    df = dfg[fields_to_aggregate].median()
    fig = px.scatter(
        df.reset_index(),
        x="log2_fold_change",
        y="-log10_pval_adjusted_bh",
        facet_col="log2_fc",
        facet_row="malignant_means",
        color="gene_perturbed",
        symbol="gene_perturbed",
        symbol_map={False: "circle", True: "x"},
    )
    fig.update_traces(marker_size=1)
    fig.update_xaxes(range=[-2, 2])
    fig.update_yaxes(range=[0, 6])
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    fig = _util_remove_excess_axis_titles(fig)
    fig.update_layout(width=1000, height=1000)
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


def calculate_roc(df: pd.DataFrame, score_column: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    # df = df.loc[df.index.get_level_values("log2_fc") != "0.00"]
    # df = df.loc[df.index.get_level_values("run_id").isin(["00", "01"])]
    every_index_level_except_gene_symbol = list(
        filter(lambda x: x != "gene_symbol", df.index.names)
    )
    logger.debug("grouping by %s", every_index_level_except_gene_symbol)
    dfg = df.groupby(every_index_level_except_gene_symbol)

    def compute_roc_curve_for_group(df: pd.DataFrame) -> pd.DataFrame:
        fpr, tpr, scores = sklearn.metrics.roc_curve(
            y_true=df["gene_perturbed"],
            y_score=df[score_column],
        )
        return pd.DataFrame({"fpr": fpr, "tpr": tpr, score_column: scores})

    def compute_roc_auc_score_for_group(df: pd.DataFrame) -> float:
        try:
            return sklearn.metrics.roc_auc_score(
                y_true=df["gene_perturbed"],
                y_score=df[score_column],
            )
        except ValueError:
            return np.nan

    roc_curves = dfg.apply(compute_roc_curve_for_group)
    roc_auc_scores = dfg.apply(compute_roc_auc_score_for_group)
    return roc_curves, roc_auc_scores


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


def calculate_precision_and_recall(df: pd.DataFrame, score_column: str) -> pd.DataFrame:
    # df = df.loc[df.index.get_level_values("log2_fc") != "0.00"]
    every_index_level_except_gene_symbol = list(
        filter(lambda x: x != "gene_symbol", df.index.names)
    )
    logger.debug("grouping by %s", every_index_level_except_gene_symbol)
    dfg = df.groupby(every_index_level_except_gene_symbol)

    def compute_precision_recall_curve_for_group(df: pd.DataFrame) -> pd.DataFrame:
        precision, recall, scores = sklearn.metrics.precision_recall_curve(
            y_true=df["gene_perturbed"],
            probas_pred=df[score_column],
        )
        # extend scores by one to include infinity
        scores = np.append(scores, np.inf)
        return pd.DataFrame(
            {
                "precision": precision,
                "recall": recall,
                score_column: scores,
            }
        )

    def compute_precision_for_group(df: pd.DataFrame) -> pd.DataFrame:
        return sklearn.metrics.precision_score(
            y_true=df["gene_perturbed"],
            y_pred=df["significant_bh_fdr=0.10"],
        )

    df_curves = dfg.apply(compute_precision_recall_curve_for_group)
    df_precision = dfg.apply(compute_precision_for_group)
    return df_curves, df_precision


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


def compute_pval_adjusted_for_group(df: pd.DataFrame) -> pd.DataFrame:
    alpha = 0.1  # false discovery rate
    significance_column = f"significant_bh_fdr={alpha:.2f}"
    reject, pvals_corrected = multipletests(df["pval"], alpha=alpha, method="fdr_bh")[0:2]
    return pd.DataFrame(
        {
            significance_column: reject,
            "pval_adjusted_bh": pvals_corrected,
            "-log10_pval_adjusted_bh": -np.log10(pvals_corrected),
        },
        index=df.index,
    )


def compute_pval_adjusted_fields(df_gene_stats: pd.DataFrame) -> pd.DataFrame:
    observed_log2_fc = df_gene_stats["log2_fold_change"].fillna(0)
    sign_observed_log2_fc = np.sign(observed_log2_fc).replace({0: 1})
    dist_log2_fc = df_gene_stats.index.to_frame()["log2_fc"].astype(float)
    sign_dist_log2_fc = np.sign(dist_log2_fc).replace({0: 1})
    # get list of every index level name except "gene_symbol"
    groupby_fields = [x for x in df_gene_stats.index.names if x != "gene_symbol"]
    logger.debug("Grouping by %s", groupby_fields)
    # don't prepend groupby fields to the resulting index
    dfg = df_gene_stats.groupby(groupby_fields, group_keys=False)
    # sanity check - assert size of each group is 5000
    logger.debug("Each group's size: %s", dfg.size())
    # assert (dfg.size() == 16063).all(), "Each group should have 16063 genes"
    df_pval_adjusted = dfg.apply(compute_pval_adjusted_for_group)
    assert df_gene_stats.index.equals(df_pval_adjusted.index), "Index of result is different"
    df_pval_adjusted["pval_adjusted_bh_signed"] = (
        df_pval_adjusted["pval_adjusted_bh"] * sign_observed_log2_fc
    )
    df_pval_adjusted["pval_adjusted_bh_signed_directional"] = (
        df_pval_adjusted["pval_adjusted_bh"] * sign_dist_log2_fc * sign_observed_log2_fc
    )
    df_pval_adjusted["-log10_pval_adjusted_bh_signed"] = (
        df_pval_adjusted["-log10_pval_adjusted_bh"] * sign_observed_log2_fc
    )
    df_pval_adjusted["-log10_pval_adjusted_bh_signed_directional"] = (
        df_pval_adjusted["-log10_pval_adjusted_bh"] * sign_dist_log2_fc * sign_observed_log2_fc
    )
    return df_pval_adjusted


def add_pval_adjusted_fields(
    df_gene_stats: pd.DataFrame, df_pval_adjusted: pd.DataFrame
) -> pd.DataFrame:
    assert df_gene_stats.index.equals(df_pval_adjusted.index), "Index mismatch"
    # fields in df_pval_adjusted that are not in df_gene_stats
    fields = filter(lambda x: x not in df_gene_stats.columns, df_pval_adjusted.columns)
    additional_fields = df_pval_adjusted[fields]
    logger.debug("adding fields %s", list(additional_fields.columns))
    return df_gene_stats.join(additional_fields, how="inner", validate="1:1")


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


def _util_remove_excess_axis_titles(fig: go.Figure) -> go.Figure:
    # remove axis titles from all subplots... except the middle one
    import statistics

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
