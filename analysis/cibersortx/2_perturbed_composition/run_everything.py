import logging
import itertools

import dask.dataframe as dd
import numpy as np
import pandas as pd
from upath import UPath

import helpers
from helpers.simulating_bulk_rnaseq import perturb_malignant_fractions
from helpers import datasets
from helpers.simulating_bulk_rnaseq.stuff_for_everything import simulate_data
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


if __name__ == "__main__":
    helpers.logging.configure_logging()
    logger.setLevel("DEBUG")

    n = 50  # samples per group
    root_all_results = UPath(
        f"gs://liulab/differential_composition/{make_a_nice_timestamp_of_now()}"
    )
    seed_counter = itertools.count()
    malignant_fraction_mean_pairs = [
        (0.71, 0.71),
        (0.7, 0.72),
        (0.65, 0.75),
        (0.6, 0.8),
        (0.55, 0.85),
    ]

    logger.debug("reading source data")
    sc_rnaseq, sc_metadata = datasets.jerby_arnon.load_scrnaseq_and_filter_genes()
    # sc_rnaseq = sc_rnaseq.iloc[::1000, :]
    f_tcga_skcm_mets = datasets.tcga_skcm.load_fractions_mets_only()
    logger.debug("shapes: %s, %s, %s", sc_rnaseq.shape, sc_metadata.shape, f_tcga_skcm_mets.shape)

    for run_id, malignant_fraction_mean_pair in itertools.product(
        range(10), malignant_fraction_mean_pairs
    ):
        malignant_means_str = ",".join(map(str, malignant_fraction_mean_pair))
        experiment_path = (
            root_all_results / f"run_id={run_id:02d}" / f"malignant_means={malignant_means_str}"
        )
        logger.debug("experiment_path: %s", experiment_path)
        rng = np.random.default_rng(seed=next(seed_counter))
        fractions_list, bulk_rnaseq_list, cell_type_geps_list = [], [], []
        for name, mean_malignant_value in zip(["low", "high"], malignant_fraction_mean_pair):
            f_tcga_skcm_mets_perturbed = perturb_malignant_fractions(
                f_tcga_skcm_mets, mean_malignant_value
            )
            fractions, bulk_rnaseq, cell_type_geps = simulate_data(
                sc_rnaseq, sc_metadata, f_tcga_skcm_mets_perturbed, rng, n, name
            )
            logger.debug("saving data to %s", experiment_path)
            cell_type_geps.to_parquet(experiment_path / name / "cell_type_geps.parquet")
            bulk_rnaseq.to_parquet(experiment_path / name / "bulk_rnaseq.parquet")
            fractions.to_parquet(experiment_path / name / "fractions.parquet")
            fractions_list.append(fractions)
            bulk_rnaseq_list.append(bulk_rnaseq)
            cell_type_geps_list.append(cell_type_geps)

        logger.debug("running cibersortx")
        bulk_rnaseq_all, fractions_all = pd.concat(bulk_rnaseq_list, axis=1), pd.concat(
            fractions_list, axis=0
        )
        logger.debug("running cibersortx; saving to %s", experiment_path / "cibersortx")
        helpers.running_cibersortx.hires_only.run_and_upload_from_dataframes(
            bulk_rnaseq_all,
            fractions_all,
            experiment_path / "cibersortx",
        )

        logger.debug("computing stats for bulk rna-seq")
        logger.debug("computing stats and saving to %s", experiment_path / "deg_analysis")
        gene_stats_bulk = compute_stats(bulk_rnaseq_all, "low", "high")
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
        gene_stats_malignant_cibersortx = compute_stats(rnaseq_malignant_cibersortx, "low", "high")
        gene_stats_malignant_cibersortx["perturbed"] = False
        gene_stats_malignant_cibersortx.to_parquet(
            experiment_path / "deg_analysis" / "gene_stats_malignant_cibersortx.parquet"
        )
