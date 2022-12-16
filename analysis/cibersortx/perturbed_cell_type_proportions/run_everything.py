import logging

import dask.dataframe as dd
import numpy as np
import pandas as pd
from upath import UPath

import helpers
from helpers import datasets
from helpers.simulating_bulk_rnaseq import simulate_data
from helpers.useful_small_things import make_a_nice_timestamp_of_now

logger = logging.getLogger(__name__)


def compute_stats(
    rnaseq: pd.DataFrame,
    group_1: str,
    group_2: str,
) -> pd.DataFrame:
    rnaseq = rnaseq.copy()
    rnaseq.columns = pd.MultiIndex.from_tuples(
        rnaseq.columns.str.split("/", expand=True), names=["group_id", "sample_id"]
    )
    rnaseq = rnaseq.stack(["group_id", "sample_id"])
    rnaseq_groupby = rnaseq.groupby("gene_symbol")
    gene_stats = helpers.deg_analysis.compute_stats(rnaseq_groupby, "group_id", group_1, group_2)
    return gene_stats


def perturb_malignant_fractions(
    fractions: pd.DataFrame, malignant_fraction_log2_fc: float
) -> pd.DataFrame:
    fractions = fractions.copy()
    factor = 2**malignant_fraction_log2_fc
    logger.debug("perturbing malignant fraction by %s", factor)
    fractions["Malignant"] *= factor
    fractions = fractions.div(fractions.sum(axis=1), axis=0)
    return fractions


if __name__ == "__main__":
    helpers.logging.configure_logging()
    logger.setLevel("DEBUG")
    logging.getLogger("pandas").setLevel("DEBUG")
    logging.getLogger("upath").setLevel("DEBUG")
    logging.getLogger("helpers").setLevel("DEBUG")
    logging.getLogger("helpers.deg_analysis").setLevel("INFO")
    logging.getLogger("helpers.simulating_bulk_rnaseq").setLevel("INFO")

    N = 50
    rng = np.random.default_rng(seed=0)
    path_results = UPath(f"gs://liulab/differential_composition/{make_a_nice_timestamp_of_now()}")
    malignant_fraction_log2_fc_values = [-2.0, -1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 2.0] * 3

    logger.debug("reading source data")
    sc_rnaseq, sc_metadata = datasets.jerby_arnon.load_scrnaseq_and_filter_genes()
    fractions_tcga_skcm_mets = datasets.tcga_skcm.load_fractions_mets_only()
    logger.debug("determining genes to perturb")

    # control data
    logger.debug("creating control data")
    fractions_control, bulk_rnaseq_control, cell_type_geps_control = simulate_data(
        sc_rnaseq, sc_metadata, fractions_tcga_skcm_mets, rng, N, name="control"
    )
    logger.debug("saving data to %s", path_results / "control")
    cell_type_geps_control.to_parquet(path_results / "control" / "cell_type_geps.parquet")
    bulk_rnaseq_control.to_parquet(path_results / "control" / "bulk_rnaseq.parquet")
    fractions_control.to_parquet(path_results / "control" / "fractions.parquet")

    # experiments
    for i, malignant_fraction_log2_fc in enumerate(malignant_fraction_log2_fc_values):
        rng_seed = i + 1
        experiment_name = (
            f"rng_seed={rng_seed:02d},malignant_frac_log2={malignant_fraction_log2_fc:.3f}"
        )
        rng = np.random.default_rng(seed=rng_seed)
        logger.debug("starting experiment %s", experiment_name)
        experiment_path = path_results / experiment_name
        assert np.allclose(
            fractions_tcga_skcm_mets,
            perturb_malignant_fractions(fractions_tcga_skcm_mets, 0.0),
        )
        fractions_perturbed = perturb_malignant_fractions(
            fractions_tcga_skcm_mets, malignant_fraction_log2_fc
        )
        logger.debug("simulating data")
        fractions_test, bulk_rnaseq_test, cell_type_geps = simulate_data(
            sc_rnaseq, sc_metadata, fractions_perturbed, rng, N, name=experiment_name
        )
        logger.debug("saving data to %s", experiment_path)
        cell_type_geps.to_parquet(experiment_path / "cell_type_geps.parquet")
        bulk_rnaseq_test.to_parquet(experiment_path / "bulk_rnaseq.parquet")
        fractions_test.to_parquet(experiment_path / "fractions.parquet")
        logger.debug("running cibersortx")
        bulk_rnaseq_all = pd.concat([bulk_rnaseq_control, bulk_rnaseq_test], axis=1)
        fractions_all = pd.concat([fractions_control, fractions_test], axis=0)
        logger.debug("running cibersortx and saving results to %s", experiment_path / "cibersortx")
        helpers.running_cibersortx.hires_only.run_and_upload_from_dataframes(
            bulk_rnaseq_all,
            fractions_all,
            experiment_path / "cibersortx",
        )

        logger.debug("computing stats for bulk rna-seq")
        logger.debug("computing stats and saving to %s", experiment_path / "deg_analysis")
        gene_stats_bulk = compute_stats(bulk_rnaseq_all, "control", experiment_name)
        gene_stats_bulk["perturbed"] = False
        gene_stats_bulk.to_parquet(experiment_path / "deg_analysis" / "gene_stats_bulk.parquet")

        logger.debug("computing stats for malignant rna-seq inferred by CIBERSORTx")
        pattern = experiment_path / "cibersortx" / "**" / "*Malignant*txt"
        logger.debug("reading cibersortx inferred rnaseq for malignant cells from %s", pattern)
        rnaseq_malignant_cibersortx = (
            dd.read_csv(pattern, sep="\t")
            .rename(columns={"GeneSymbol": "gene_symbol"})
            .set_index("gene_symbol")
            .compute()
        )
        gene_stats_malignant_cibersortx = compute_stats(
            rnaseq_malignant_cibersortx, "control", experiment_name
        )
        gene_stats_malignant_cibersortx["perturbed"] = False
        gene_stats_malignant_cibersortx.to_parquet(
            experiment_path / "deg_analysis" / "gene_stats_malignant_cibersortx.parquet"
        )
