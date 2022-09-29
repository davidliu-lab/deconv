import logging
import warnings

import numpy as np
import pandas as pd
import plotly.basedatatypes
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import scipy.stats
from pandas.core.groupby.generic import SeriesGroupBy
from statsmodels.stats.multitest import multipletests

logger = logging.getLogger(__name__)


def calculate_pval_threshold(
    test_was_fdr_significant: pd.Series, alpha: float
) -> float:
    n_signif_results = test_was_fdr_significant.sum()
    pval_threshold = (n_signif_results + 1) * alpha / len(test_was_fdr_significant)
    return pval_threshold


def add_multipletests_stats(df: pd.DataFrame) -> pd.DataFrame:
    alphas = [0.1, 0.25]  # false discovery rates
    df["-log10_pval"] = -np.log10(df["pval"])
    df["log2_fold_change"] = np.log2(df["fold_change"])
    sign = np.sign(df["log2_fold_change"])
    df["-log10_pval_signed"] = df["-log10_pval"] * sign
    # multiple hypothesis testing with benjamini-hochberg
    for alpha in alphas:
        significance_column = f"significant_bh_fdr={alpha:.2f}"
        df[significance_column] = multipletests(
            df["pval"], alpha=alpha, method="fdr_bh"
        )[0]
        n_signif_results = df[significance_column].sum()
        pval_threshold_str = f"pval_threshold_fdr={alpha:.2f}"
        df.attrs[pval_threshold_str] = (n_signif_results + 1) * alpha / len(df)
        df.attrs[f"-log10_{pval_threshold_str}"] = -np.log10(
            df.attrs[pval_threshold_str]
        )
    return df


def compute_stats(
    series_groupby: SeriesGroupBy,
    group_col: str,
    group_1: str,
    group_2: str,
) -> pd.DataFrame:
    def compute_stats_for_group(series: pd.Series) -> pd.Series:
        logger.debug("computing stats for group %s", series.name)
        series_groupby = series.groupby(group_col)
        series_1 = series_groupby.get_group(group_1)
        series_2 = series_groupby.get_group(group_2)
        logger.debug("shapes: %s", (series_1.shape, series_2.shape))
        sparsity_overall = np.mean(series == 0)
        logger.debug("sparsity overall: %s", sparsity_overall)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # ignore divide by zero warnings
            fold_change = series_2.mean() / series_1.mean()
        results = pd.Series(
            {
                "pval": scipy.stats.mannwhitneyu(series_1, series_2)[1],
                # "pval_ttest": scipy.stats.ttest_ind(series_1, series_2).pvalue,
                "fold_change": fold_change,
                "sparsity_overall": sparsity_overall,
            }
        )
        return results

    df_stats_by_group = series_groupby.apply(compute_stats_for_group)
    df_stats_by_group = df_stats_by_group.unstack(-1)
    df_stats_by_group = df_stats_by_group.reset_index()
    df_stats_by_group = add_multipletests_stats(df_stats_by_group)
    return df_stats_by_group


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
            annotation_text=f"FDR = {alpha:.2f}",
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
            annotation_text=f"FDR = {alpha:.2f}",
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


def make_volcano_facets(path_root, log2_scalars, horizontal=True):
    n_scalars = len(log2_scalars)
    files = dict()
    for log2_scaling_factor in log2_scalars:
        scaling_factor = 2**log2_scaling_factor
        scaling_factor_str = f"scaling_factor={scaling_factor:.3f}"
        files[(scaling_factor, "Bulk simulated")] = pd.read_parquet(
            path_root / scaling_factor_str / "gene_stats_bulk.parquet"
        )
        files[(scaling_factor, "Malignant inferred")] = pd.read_parquet(
            path_root / scaling_factor_str / "gene_stats_malignant.parquet"
        )
    df = pd.concat(files, names=["scaling_factor", "data_origin"]).reset_index()
    df = df.rename(
        columns={
            "log2_fold_change": r"$\log_{2} [\text{fold change}]$",
            "-log10_pval": r"$-\log_{10} [\text{p-value}]$",
        }
    )
    fig = px.scatter(
        df,
        x=r"$\log_{2} [\text{fold change}]$",
        y=r"$-\log_{10} [\text{p-value}]$",
        color="perturbed",
        # color=df["perturbed"].map({True: "red", False: "blue"}),
        # size=df["perturbed"].map({True: 2, False: 1.5}),
        facet_col="scaling_factor" if horizontal else "data_origin",
        facet_row="data_origin" if horizontal else "scaling_factor",
        hover_name="gene_symbol",
        hover_data=[
            "gene_symbol",
            "perturbed",
            "pval",
            "fold_change",
            "sparsity_overall",
        ],
    )
    fig.update_xaxes(range=[-8, 8])
    fig.update_yaxes(range=[0, 18])
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("data_origin=")[-1]))
    for k, log2_scaling_factor in enumerate(log2_scalars):
        for i, data_origin in enumerate(["Bulk simulated", "Malignant inferred"]):
            row_number = (
                i + 1 if horizontal else n_scalars - k
            )  # bug with plotly row selection!
            column_number = k + 1 if horizontal else i + 1
            logger.debug(
                "adding data for k=%d, row_number=%d, column_number=%d",
                k,
                row_number,
                column_number,
            )
            scaling_factor = 2**log2_scaling_factor
            for alpha in [0.1, 0.25]:
                key = (scaling_factor, data_origin)
                significances = files[key][f"significant_bh_fdr={alpha:.2f}"]
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
                    annotation_text=f"FDR = {alpha:.2f}",
                    annotation_position="top left" if alpha == 0.1 else "bottom left",
                    row=row_number,
                    col=column_number,
                )
    fig.update_layout(
        height=800 if horizontal else 300 * n_scalars,
        width=n_scalars * 300 if horizontal else 1000,
        font=dict(family="Courier New, monospace"),
        legend=dict(yanchor="top", y=0.95, xanchor="right", x=0.95),
        # marker=dict(line=dict(width=0)),
    )
    fig.update_xaxes(dtick=2.0)  # , col="all")
    fig.update_xaxes(showticklabels=True)
    fig.update_yaxes(showticklabels=True)
    return fig
