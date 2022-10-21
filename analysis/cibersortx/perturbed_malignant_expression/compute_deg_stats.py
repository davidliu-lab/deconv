import logging
from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd
from upath import UPath

import helpers

logger = logging.getLogger(__name__)


def load_data(path: Union[UPath, Path]) -> pd.Series:
    rows_to_skip = []
    # rows_to_skip = lambda i: i % 25
    df = pd.read_csv(path, sep="\t", index_col=0, skiprows=rows_to_skip)
    df = df.rename_axis(index="gene_symbol")
    df.columns = pd.MultiIndex.from_tuples(
        df.columns.str.split("/", expand=True), names=["group_id", "sample_id"]
    )
    return df.stack(["group_id", "sample_id"])


def compute_for(path: Union[UPath, Path], scaling_factor_str: str) -> pd.DataFrame:
    logger.debug("reading %s", path)
    series_rnaseq = load_data(path)
    groups = series_rnaseq.groupby("gene_symbol")
    df_gene_stats = helpers.deg_analysis.compute_stats(
        groups, "group_id", "control", scaling_factor_str
    )
    return df_gene_stats


if __name__ == "__main__":
    timestamp_str = helpers.useful_small_things.make_a_nice_timestamp_of_now()
    SAVE = True
    logger.setLevel("DEBUG")
    logging.getLogger("helpers").setLevel("DEBUG")
    logging.getLogger("helpers.deg_analysis").setLevel("INFO")
    logging.basicConfig(format=helpers.logging.format_string)
    genes_perturbed = pd.read_csv(
        "gs://liulab/simulated/perturbed_malignant_expression/20221015_21h52m40s/genes_perturbed.csv"
    )["gene_symbol"]
    path_root_cibersortx_results = (
        UPath("gs://liulab/cibersortx/perturbed_malignant_expression")
        / "20221019_15h23m14s"
    )

    paths_cibersortx_results = [
        p.parent for p in path_root_cibersortx_results.glob("**/outdir")
    ]
    logger.debug("paths: %s", paths_cibersortx_results)
    path_target_root = UPath("gs://liulab/deg_analysis") / timestamp_str
    logger.debug("path_target_root: %s", path_target_root)
    for path_cibersortx_result in paths_cibersortx_results:
        logger.debug("processing %s", path_cibersortx_result)
        result_description = path_cibersortx_result.name
        assert result_description.startswith("log2_fc="), result_description
        # bulk simulated
        logger.debug("computing for bulk")
        path_bulk_rnaseq = next(
            (path_root_cibersortx_results / result_description).glob(
                "**/bulkrnaseq.txt"
            )
        )
        df_gene_stats_bulk = compute_for(path_bulk_rnaseq, result_description)
        df_gene_stats_bulk["perturbed"] = df_gene_stats_bulk["gene_symbol"].isin(
            genes_perturbed
        )
        # logger.debug("df_gene_stats_bulk: %s", df_gene_stats_bulk)
        if SAVE:
            path_target_bulk = (
                path_target_root / result_description / "gene_stats_bulk.parquet"
            )
            logger.debug("writing %s", path_target_bulk)
            df_gene_stats_bulk.to_parquet(path_target_bulk)

        # malignant inferred
        logger.debug("computing for inferred malignant")
        path_malignant = next(
            (path_root_cibersortx_results / result_description).glob(
                "**/CIBERSORTxHiRes_NA_Malignant*txt"
            )
        )
        df_gene_stats_malignant = compute_for(path_malignant, result_description)
        df_gene_stats_malignant["perturbed"] = df_gene_stats_malignant[
            "gene_symbol"
        ].isin(genes_perturbed)
        # logger.debug("df_gene_stats_malignant: %s", df_gene_stats_malignant)
        if SAVE:
            path_target_malignant = (
                path_target_root / result_description / "gene_stats_malignant.parquet"
            )
            logger.debug("writing %s", path_target_malignant)
            df_gene_stats_malignant.to_parquet(path_target_malignant)
