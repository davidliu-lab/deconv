import itertools
import logging
from typing import Union

import dask.dataframe as dd
import numpy as np
import pandas as pd
from upath import UPath

import helpers
from helpers import datasets
from helpers.simulating_bulk_rnaseq import perturb_malignant_fractions
from helpers.simulating_bulk_rnaseq.gene_perturbation import (
    select_100_genes_at_least_somewhat_expressed_in_malignant,
    perturb_scrnaseq_gene_expression,
)
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
        rnaseq.columns.str.split("/", expand=True),
        names=["group_id", "sample_id"],
    )
    rnaseq = rnaseq.stack(["group_id", "sample_id"])
    rnaseq_groupby = rnaseq.groupby("gene_symbol")
    gene_stats = helpers.deg_analysis.compute_stats(rnaseq_groupby, "group_id", group_1, group_2)
    return gene_stats


def construct_dataset(
    sc_rnaseq: pd.DataFrame,
    sc_metadata: pd.DataFrame,
    fractions_to_sample_from: pd.DataFrame,
    genes_to_perturb: pd.Index,
    log2_fc: Union[float, None],
    mean_malignant_value: Union[float, None],
    rng: np.random.Generator,
    n: int,
    name: str,
):
    if log2_fc:
        sc_rnaseq = perturb_scrnaseq_gene_expression(
            sc_rnaseq,
            sc_metadata,
            "Malignant",
            genes_to_perturb,
            scaling_factor=2.0**log2_fc,
        )
    if mean_malignant_value:
        fractions_to_sample_from = perturb_malignant_fractions(
            fractions_to_sample_from, mean_malignant_value
        )
    else:
        fractions_to_sample_from = fractions_to_sample_from
    fractions, bulk_rnaseq, cell_type_geps = simulate_data(
        sc_rnaseq,
        sc_metadata,
        fractions_to_sample_from,
        rng,
        n,
        name,
    )
    return fractions, bulk_rnaseq, cell_type_geps


if __name__ == "__main__":
    helpers.logging.configure_logging()
    logger.setLevel("DEBUG")
    n = 50  # samples per group
    root_results = UPath(
        f"gs://liulab/differential_composition_and_expression/{make_a_nice_timestamp_of_now()}"
    )
    seed_counter = itertools.count()
    experiment_counter = itertools.count()
    rng = np.random.default_rng(seed=next(seed_counter))

    logger.debug("reading source data")
    (
        sc_rnaseq,
        sc_metadata,
    ) = datasets.jerby_arnon.load_scrnaseq_and_filter_genes()
    # sc_rnaseq = sc_rnaseq.iloc[::100, :]  # faster debug
    f_tcga_skcm_mets = datasets.tcga_skcm.load_fractions_mets_only()
    logger.debug(
        "shapes: %s, %s, %s",
        sc_rnaseq.shape,
        sc_metadata.shape,
        f_tcga_skcm_mets.shape,
    )
    genes_to_perturb = select_100_genes_at_least_somewhat_expressed_in_malignant(
        sc_rnaseq, sc_metadata, rng
    )
    malignant_log2_fold_changes = np.array([-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5])
    malignant_fraction_mean_pairs = [
        (0.55, 0.85),
        (0.65, 0.75),
        (0.7, 0.72),
        (None, None),
        (0.71, 0.71),
        (0.72, 0.7),
        (0.75, 0.65),
        (0.85, 0.55),
    ]
    pd.DataFrame(genes_to_perturb).to_csv(root_results / "genes_perturbed.csv", index=False)
    for run_id, malignant_fraction_mean_pair, log2_fc in itertools.product(
        range(5), malignant_fraction_mean_pairs, malignant_log2_fold_changes
    ):
        logger.debug("starting experiment")
        malignant_means_str = ",".join(map(str, malignant_fraction_mean_pair))
        experiment_path = (
            root_results
            / f"experiment_id={next(experiment_counter):03d}"
            / f"malignant_means={malignant_means_str}"
            / f"log2_fc={log2_fc:4.2f}"
            / f"run_id={run_id:02d}"
        )
        logger.debug("experiment_path: %s", experiment_path)
        rng = np.random.default_rng(seed=next(seed_counter))
        logger.debug("perturbing gene expression")
        fractions_list, bulk_rnaseq_list = [], []
        for name, mean_malignant_value in zip(["a", "b"], malignant_fraction_mean_pair):
            log2_fc = log2_fc if name == "b" else 0.0
            fractions, bulk_rnaseq, cell_type_geps = construct_dataset(
                sc_rnaseq,
                sc_metadata,
                f_tcga_skcm_mets,
                genes_to_perturb,
                log2_fc,
                mean_malignant_value,
                rng,
                n,
                name,
            )
            logger.debug("saving data to %s", experiment_path)
            cell_type_geps.to_parquet(experiment_path / name / "cell_type_geps.parquet")
            bulk_rnaseq.to_parquet(experiment_path / name / "bulk_rnaseq.parquet")
            fractions.to_parquet(experiment_path / name / "fractions.parquet")
            fractions_list.append(fractions)
            bulk_rnaseq_list.append(bulk_rnaseq)

        logger.debug("running cibersortx")
        bulk_rnaseq_all = pd.concat(bulk_rnaseq_list, axis=1)
        fractions_all = pd.concat(fractions_list, axis=0)
        logger.debug("running cibersortx; saving to %s", experiment_path / "cibersortx")
        helpers.running_cibersortx.hires_only.run_and_upload_from_dataframes(
            bulk_rnaseq_all,
            fractions_all,
            experiment_path / "cibersortx",
        )

        logger.debug("computing stats for bulk rna-seq")
        logger.debug("computing stats and saving to %s", experiment_path / "deg_analysis")
        gene_stats_bulk = compute_stats(bulk_rnaseq_all, "a", "b")
        gene_stats_bulk["perturbed"] = gene_stats_bulk["gene_symbol"].isin(genes_to_perturb)
        gene_stats_bulk.to_parquet(experiment_path / "deg_analysis" / "gene_stats_bulk.parquet")

        logger.debug("computing stats for malignant rna-seq inferred by CIBERSORTx")
        pattern = experiment_path / "cibersortx" / "**" / "*Malignant*txt"
        logger.debug(
            "reading cibersortx inferred rnaseq for malignant cells from %s",
            pattern,
        )
        rnaseq_malignant_cibersortx = (
            dd.read_csv(pattern, sep="\t")
            .rename(columns={"GeneSymbol": "gene_symbol"})
            .set_index("gene_symbol")
            .compute()
        )
        gene_stats_malignant_cibersortx = compute_stats(rnaseq_malignant_cibersortx, "a", "b")
        gene_stats_malignant_cibersortx["perturbed"] = False
        gene_stats_malignant_cibersortx.to_parquet(
            experiment_path / "deg_analysis" / "gene_stats_malignant_cibersortx.parquet"
        )
