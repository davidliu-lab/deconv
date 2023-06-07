import logging

import numpy as np
import pandas as pd
import plotly.basedatatypes
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from regex import F

from helpers.deg_analysis import plotting_utils
from helpers.deg_analysis.pvals_with_fdr import calculate_pval_threshold
from helpers.deg_analysis.plotting_utils import _util_remove_excess_axis_titles

logger = logging.getLogger(__name__)


def make_volcano_figure(df_stats: pd.DataFrame, perturbed: bool = False):
    fig = px.scatter(
        df_stats,
        x="log2_fold_change",
        y="-log10_pval",
        color="perturbed" if perturbed else None,
        hover_name="gene_symbol",
        hover_data=["pval", "sparsity_overall"],
    )
    fig.update_layout(
        xaxis_title=r"$\log_{2} [\text{fold change}]$",
        yaxis_title=r"$-\log_{10} [\text{p-value (Mann-Whitney U)}]$",
        legend_title="Perturbed?" if perturbed else None,
        font=dict(family="Courier New, monospace", color="RebeccaPurple"),
        height=500,
        width=500,
    )
    fig.update_traces(marker=dict(size=3))
    fig.update_layout(legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01))
    for alpha in [0.1, 0.25]:
        significances = df_stats[f"significant_bh_fdr={alpha:.2f}"]
        pval_threshold = calculate_pval_threshold(significances, alpha)
        fig.add_hline(
            y=-np.log10(pval_threshold),
            line_dash="dot",
            annotation_text=f"FDR={alpha:.2f}",
            annotation_position="top left",
        )
    return fig


def make_volcano_trace(df_stats: pd.DataFrame) -> plotly.basedatatypes.BaseTraceType:
    trace = go.Scattergl(
        x=df_stats["log2_fold_change"],
        y=df_stats["-log10_pval"],
        mode="markers",
        marker=dict(
            size=10,
            color=df_stats["perturbed"].astype(int),
            colorscale=["black", "blue"],
            opacity=0.5,
            symbol=df_stats["perturbed"].map({True: "star", False: "circle"}),
        ),
        hovertext=df_stats["gene_symbol"],
        customdata=df_stats[["pval", "sparsity_overall"]],
        hovertemplate="Gene: %{hovertext}<br>p-value: %{customdata[0]:.2e}<br>Sparsity: %{customdata[1]:.2f}<extra></extra>",
    )
    return trace


def make_volcano_traces(
    df_stats: pd.DataFrame,
) -> list[plotly.basedatatypes.BaseTraceType]:
    traces = []
    df_groupby = df_stats.groupby("perturbed")
    for perturbed, df_group in df_groupby:
        logger.debug("making volcano trace for perturbed=%s", perturbed)
        marker_kwargs = dict(
            color="red" if perturbed else "blue",
            size=6 if perturbed else 3,
            # symbol="asterisk" if perturbed else "circle",
        )
        traces.append(
            go.Scattergl(
                x=df_group["log2_fold_change"],
                y=df_group["-log10_pval"],
                mode="markers",
                marker=marker_kwargs,
                hovertext=df_stats["gene_symbol"],
                customdata=df_stats[["pval", "sparsity_overall"]],
                hovertemplate=(
                    "Gene: %{hovertext}"
                    "<br>p-value: %{customdata[0]:.2e}"
                    "<br>Sparsity: %{customdata[1]:.2f}"
                    "<extra></extra>"
                ),
            )
        )
    return traces


def make_volcano_figure_2(df_stats: pd.DataFrame) -> go.Figure:
    fig = go.Figure()
    for trace in make_volcano_traces(df_stats):
        fig.add_trace(trace)
    fig.update_layout(
        xaxis_title=r"$\log_{2} [\text{fold change}]$",
        yaxis_title=r"$-\log_{10} [\text{p-value (Mann-Whitney U)}]$",
        font=dict(family="Courier New, monospace", color="RebeccaPurple"),
        height=500,
        width=500,
        showlegend=False,
    )
    fig.update_xaxes(range=[-5, 5])
    fig.update_yaxes(range=[0, 20])
    for alpha in [0.1]:
        significances = df_stats[f"significant_bh_fdr={alpha:.2f}"]
        pval_threshold = calculate_pval_threshold(significances, alpha)
        fig.add_hline(
            y=-np.log10(pval_threshold),
            line_dash="dot",
            annotation_text=f"FDR={alpha:.2f}",
            annotation_position="top left",
        )
    return fig


def make_volcano_subplots_figure(path_root, log2_scalars):
    n_scalars = len(log2_scalars)
    fig = make_subplots(
        cols=n_scalars,
        rows=2,
        column_titles=list(map(lambda x: f"{2**x:.3f}", log2_scalars)),
        row_titles=["Simulated bulk RNA-seq", "Inferred malignant (cibersortx)"],
    )
    for i, log2_scaling_factor in enumerate(log2_scalars):
        scaling_factor = 2**log2_scaling_factor
        scaling_factor_str = f"scaling_factor={scaling_factor:.3f}"
        # bulk simulated
        logger.debug("processing bulk simulated, %s", scaling_factor_str)
        df_gene_stats_bulk = pd.read_parquet(
            path_root / scaling_factor_str / "gene_stats_bulk.parquet"
        )
        for trace_bulk in make_volcano_traces(df_gene_stats_bulk):
            fig.add_trace(trace_bulk, row=1, col=i + 1)
        fig.update_xaxes(
            range=[-6, 6],
            title_text=r"$\log_{2} [\text{fold change}]$",
            row=1,
            col=i + 1,
        )
        fig.update_yaxes(
            range=[0, 18],
            title_text=r"$-\log_{10} [\text{p-value}]$",
            row=1,
            col=i + 1,
        )
        # malignant from cibersortx
        logger.debug("processing malignant inferred, %s", scaling_factor_str)
        df_malignant_gene_stats = pd.read_parquet(
            path_root / scaling_factor_str / "gene_stats_malignant.parquet"
        )
        for trace_malignant in make_volcano_traces(df_malignant_gene_stats):
            fig.add_trace(trace_malignant, row=2, col=i + 1)
        fig.update_xaxes(
            range=[-6, 6],
            title_text=r"$\log_{2} [\text{fold change}]$",
            row=2,
            col=i + 1,
        )
        fig.update_yaxes(
            range=[0, 18],
            title_text=r"$-\log_{10} [\text{p-value}]$",
            row=2,
            col=i + 1,
        )
    fig.update_layout(showlegend=False, width=300 * n_scalars, height=800)
    return fig


def make_volcano_subplots_figure_vertical(path_root, log2_scalars):
    n_scalars = len(log2_scalars)
    fig = make_subplots(
        cols=2,
        rows=n_scalars,
        column_titles=["Simulated bulk RNA-seq", "Inferred malignant (cibersortx)"],
        row_titles=list(map(lambda x: f"{2**x:.3f}", log2_scalars)),
    )
    for i, log2_scaling_factor in enumerate(log2_scalars):
        scaling_factor = 2**log2_scaling_factor
        scaling_factor_str = f"scaling_factor={scaling_factor:.3f}"
        # bulk simulated
        logger.debug("processing bulk simulated, %s", scaling_factor_str)
        df_gene_stats_bulk = pd.read_parquet(
            path_root / scaling_factor_str / "gene_stats_bulk.parquet"
        )
        for trace_bulk in make_volcano_traces(df_gene_stats_bulk):
            fig.add_trace(trace_bulk, row=i + 1, col=1)
        fig.update_xaxes(
            range=[-6, 6],
            title_text=r"$\log_{2} [\text{fold change}]$",
            row=i + 1,
            col=1,
        )
        fig.update_yaxes(
            range=[0, 18],
            title_text=r"$-\log_{10} [\text{p-value (Mann-Whitney U)}]$",
            row=i + 1,
            col=1,
        )
        # malignant from cibersortx
        logger.debug("processing malignant inferred, %s", scaling_factor_str)
        df_malignant_gene_stats = pd.read_parquet(
            path_root / scaling_factor_str / "gene_stats_malignant.parquet"
        )
        for trace_malignant in make_volcano_traces(df_malignant_gene_stats):
            fig.add_trace(trace_malignant, row=i + 1, col=2)
        fig.update_xaxes(
            range=[-6, 6],
            title_text=r"$\log_{2} [\text{fold change}]$",
            row=i + 1,
            col=2,
        )
        fig.update_yaxes(
            range=[0, 18],
            title_text=r"$-\log_{10} [\text{p-value (Mann-Whitney U)}]$",
            row=i + 1,
            col=2,
        )
    fig.update_layout(showlegend=False, width=800, height=400 * n_scalars)
    return fig


def make_volcano_facets(gene_stats: pd.DataFrame, horizontal: bool = True) -> go.Figure:
    # gene_stats = gene_stats.rename(
    #     columns={
    #         "log2_fold_change": r"$\log_{2} [\text{fold change}]$",
    #         "-log10_pval": r"$-\log_{10} [\text{p-value}]$",
    #     }
    # )
    fig = px.scatter(
        gene_stats,
        # x=r"$\log_{2} [\text{fold change}]$",
        # y=r"$-\log_{10} [\text{p-value}]$",
        x="log2_fold_change",
        y="-log10_pval",
        color="perturbed",
        # color=gene_stats["perturbed"].map({True: "red", False: "blue"}),
        # size=gene_stats["perturbed"].map({True: 2, False: 1.5}),
        facet_col="experiment_name" if horizontal else "data_origin",
        facet_row="data_origin" if horizontal else "experiment_name",
        hover_name="gene_symbol",
        hover_data=[
            "gene_symbol",
            "perturbed",
            "pval",
            "fold_change",
            "sparsity_overall",
        ],
        labels={
            "log2_fold_change": r"$\log_{2} [\text{fold change}]$",
            "-log10_pval": r"$-\log_{10} [\text{p-value}]$",
        },
    )
    fig.update_xaxes(range=[-8, 8])
    fig.update_yaxes(range=[0, 18])
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("data_origin=")[-1]))
    n_experiments = len(gene_stats["experiment_name"].unique())
    fig.update_layout(
        height=800 if horizontal else 400 * n_experiments,
        width=n_experiments * 400 if horizontal else 1000,
        font=dict(family="Courier New, monospace"),
        legend=dict(yanchor="top", y=0.95, xanchor="right", x=0.95),
        # marker=dict(line=dict(width=0)),
    )
    fig.update_xaxes(dtick=2.0)  # , col="all")
    fig.update_xaxes(showticklabels=True)
    fig.update_yaxes(showticklabels=True)
    return fig


def add_fdr_lines_old(fig, gene_stats, horizontal=True, alpha_values=(0.1, 0.25)):
    experiments = gene_stats["experiment_name"].unique()
    n_experiments = len(experiments)
    for k, experiment_name in enumerate(experiments):
        for i, data_origin in enumerate(["Bulk simulated", "Malignant inferred"]):
            row_number = (
                i + 1 if horizontal else n_experiments - k
            )  # bug with plotly row selection!
            column_number = k + 1 if horizontal else i + 1
            logger.debug(
                "adding data for k=%d, row_number=%d, column_number=%d",
                k,
                row_number,
                column_number,
            )
            for alpha in alpha_values:
                key = (experiment_name, data_origin)
                experiment_gene_stats = gene_stats[gene_stats["experiment_name"] == experiment_name]
                significances = experiment_gene_stats[f"significant_bh_fdr={alpha:.2f}"]
                pval_threshold = calculate_pval_threshold(significances, alpha)
                logger.debug(
                    "-log10(pval_threshold) for alpha=%s, key=%s:, %s",
                    alpha,
                    key,
                    -np.log10(pval_threshold),
                )
                fig.add_hline(
                    y=-np.log10(pval_threshold),
                    line_dash="dot",
                    annotation_text=f"FDR={alpha:.2f}",
                    annotation_position="top left" if alpha == 0.1 else "bottom left",
                    row=row_number,
                    col=column_number,
                )
    return fig


def make_volcano_grid_scatter(
    df: pd.DataFrame,
    groupby_cols: list[str] = [],
    pval_col: str = "-log10_pval_adjusted_bh",
    hover_data: list[str] = [],
    perturbed_col: str = "gene_perturbed",
    facet_col: str = "log2_fc",
    facet_row: str = "malignant_means",
    simplify_facet_titles: bool = True,
    with_fdr_lines: bool = True,
    marker_color: str = "",
) -> go.Figure:
    fields = set(df.columns) | set(df.index.names)
    assert pval_col in fields, f"{pval_col} not in {fields}"
    assert perturbed_col in fields, f"{perturbed_col} not in {fields}"
    if groupby_cols:
        assert all([col in fields for col in groupby_cols])
        assert perturbed_col in groupby_cols, f"{perturbed_col} not in {groupby_cols}"
    else:
        if len(df.index.names) < 2:
            raise ValueError("dataframe doesn't have index levels; specify groupby_cols")
        groupby_cols = list(filter(lambda x: x != "run_id", df.index.names))
        groupby_cols.append(perturbed_col)
    logger.debug("Grouping by %s", groupby_cols)
    dfg = df.groupby(groupby_cols)
    logger.debug("Distribution of group sizes %s", dfg.size().value_counts())
    fields_to_aggregate = ["log2_fold_change", pval_col] + hover_data
    df_to_plot = dfg[fields_to_aggregate].median()
    logger.debug("plotting df_to_plot with shape %s", df_to_plot.shape)
    if not marker_color:
        marker_color = perturbed_col
    fig = px.scatter(
        df_to_plot.reset_index(),
        x="log2_fold_change",
        y=pval_col,
        facet_col=facet_col,
        facet_row=facet_row,
        color=marker_color,
        symbol=perturbed_col,
        symbol_map={False: "circle", True: "x"},
        hover_name="gene_symbol",
        hover_data=hover_data,
    )
    fig.update_traces(marker_size=5)
    fig.update_xaxes(range=[-2, 2])
    fig.update_yaxes(range=[0, 6])
    if simplify_facet_titles:
        fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
        fig = _util_remove_excess_axis_titles(fig)
    if with_fdr_lines:
        fig = plotting_utils.add_fdr_lines(fig, axis="y")
    fig.update_layout(width=800, height=500)
    # overlay legend on top right
    fig.update_layout(
        legend=dict(yanchor="top", y=0.95, xanchor="right", x=0.95),
    )
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
