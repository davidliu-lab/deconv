import itertools
import logging

import dask.dataframe as dd
import numpy as np
import pandas as pd
from upath import UPath

import helpers
from helpers import datasets
from helpers.simulating_bulk_rnaseq import perturb_malignant_fractions
from helpers.simulating_bulk_rnaseq.gene_perturbation import (
    perturb_scrnaseq_gene_expression,
    select_100_genes_at_least_somewhat_expressed_in_malignant,
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
    gene_stats = helpers.deg_analysis.stats_testing_with_fdr.compute_stats(
        rnaseq_groupby, "group_id", group_1, group_2
    )
    return gene_stats


if __name__ == "__main__":
    helpers.logging.configure_logging()
    logger.setLevel("DEBUG")
    N = 50  # samples per group
    root_results = UPath(
        f"gs://liulab/differential_composition_and_expression/{make_a_nice_timestamp_of_now()}"
    )
    seed_counter = itertools.count(start=10000)
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
    malignant_log2_fc_group_b_values = np.array([-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5])
    malignant_fraction_mean_pairs = [
        (None, None),
        (0.55, 0.85),
        (0.6, 0.8),
        (0.65, 0.75),
        (0.75, 0.65),
        (0.8, 0.6),
        (0.85, 0.55),
        (0.7, 0.9),
        (0.5, 0.7),
        (0.57, 0.83),
        (0.63, 0.77),
        (0.7, 0.72),
        (0.71, 0.71),
        (0.72, 0.7),
        (0.77, 0.63),
        (0.83, 0.57),
    ]
    pd.DataFrame(genes_to_perturb).to_csv(root_results / "genes_perturbed.csv", index=False)
    for run_id, malignant_fraction_mean_pair, log2_fc_group_b in itertools.product(
        range(20), malignant_fraction_mean_pairs, malignant_log2_fc_group_b_values
    ):
        experiment_id = next(experiment_counter)
        experiment_dict = {
            "experiment_id": experiment_id,
            # "malignant_mean_group_a": malignant_fraction_mean_pair[0],
            # "malignant_mean_group_b": malignant_fraction_mean_pair[1],
            "malignant_means": "{0},{1}".format(*malignant_fraction_mean_pair),
            "log2_fc": log2_fc_group_b,
            "run_id": run_id,
        }
        logger.debug("starting experiment")
        cibersortx_outputs_path = (
            root_results
            / "cibersortx_outputs"
            / f"experiment_id={experiment_id:03d}"
            # / f"malignant_mean_group_a={malignant_fraction_mean_pair[0]}"
            # / f"malignant_mean_group_b={malignant_fraction_mean_pair[1]}"
            / "malignant_means={0},{1}".format(*malignant_fraction_mean_pair)
            / f"log2_fc={log2_fc_group_b:4.2f}"
            / f"run_id={run_id:02d}"
        )
        # logger.debug("experiment_path: %s", experiment_path)
        logger.debug("perturbing gene expression")
        fractions_list, bulk_rnaseq_list, cell_type_geps_list = [], [], []
        for name, mean_malignant_value in zip(["a", "b"], malignant_fraction_mean_pair):
            rng = np.random.default_rng(seed=next(seed_counter))
            log2_fc = log2_fc_group_b if name == "b" else 0.0
            if log2_fc:
                logger.debug("making perturbed scRNA-seq")
                sc_rnaseq_to_sample_from = perturb_scrnaseq_gene_expression(
                    sc_rnaseq,
                    sc_metadata,
                    "Malignant",
                    genes_to_perturb,
                    scaling_factor=2.0**log2_fc,
                )
            else:
                sc_rnaseq_to_sample_from = sc_rnaseq
            if mean_malignant_value:
                logger.debug("making perturbed fractions_to_sample_from")
                fractions_to_sample_from = perturb_malignant_fractions(
                    f_tcga_skcm_mets, mean_malignant_value
                )
            else:
                fractions_to_sample_from = f_tcga_skcm_mets
            fractions, bulk_rnaseq, cell_type_geps = simulate_data(
                sc_rnaseq_to_sample_from,
                sc_metadata,
                fractions_to_sample_from,
                rng,
                N,
                name,
            )
            logger.debug("saving data")
            # cell_type_geps.to_parquet(experiment_path / name / "cell_type_geps.parquet")
            dataframes = {
                "cell_type_geps": cell_type_geps,
                "fractions": fractions,
                "bulk_rnaseq": bulk_rnaseq,
            }
            for dataset_name, df in dataframes.items():
                df = df.assign(**experiment_dict)
                df = df.assign(**{"name": name})
                partition_cols = list(experiment_dict.keys()) + ["name"]
                df.to_parquet(
                    path=root_results / dataset_name,
                    engine="pyarrow",
                    partition_cols=partition_cols,
                )
            cell_type_geps_list.append(cell_type_geps)
            fractions_list.append(fractions)
            bulk_rnaseq_list.append(bulk_rnaseq)

        logger.debug("running cibersortx")
        bulk_rnaseq_all = pd.concat(bulk_rnaseq_list, axis=1)
        fractions_all = pd.concat(fractions_list, axis=0)
        logger.debug("running cibersortx; saving to %s", cibersortx_outputs_path)
        helpers.running_cibersortx.hires_only.run_and_upload_from_dataframes(
            bulk_rnaseq_all,
            fractions_all,
            cibersortx_outputs_path,
        )

        logger.debug("computing stats for bulk rna-seq")
        # logger.debug("computing stats and saving to %s", experiment_path / "deg_analysis")
        gene_stats_bulk = compute_stats(bulk_rnaseq_all, "a", "b")
        gene_stats_bulk["perturbed"] = gene_stats_bulk["gene_symbol"].isin(genes_to_perturb)
        gene_stats_bulk = gene_stats_bulk.assign(**experiment_dict)
        gene_stats_bulk = gene_stats_bulk.assign(**{"origin": "bulk"})
        gene_stats_bulk.to_parquet(
            path=root_results / "deg_analysis",
            engine="pyarrow",
            partition_cols=list(experiment_dict.keys()) + ["origin"],
        )

        logger.debug("computing stats for malignant rna-seq inferred by CIBERSORTx")
        pattern = cibersortx_outputs_path / "**" / "*Malignant*txt"
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
        gene_stats_malignant_cibersortx["perturbed"] = gene_stats_malignant_cibersortx[
            "gene_symbol"
        ].isin(genes_to_perturb)
        gene_stats_malignant_cibersortx = gene_stats_malignant_cibersortx.assign(**experiment_dict)
        gene_stats_malignant_cibersortx = gene_stats_malignant_cibersortx.assign(
            **{"origin": "malignant_cibersortx"}
        )
        gene_stats_malignant_cibersortx.to_parquet(
            path=root_results / "deg_analysis",
            engine="pyarrow",
            partition_cols=list(experiment_dict.keys()) + ["origin"],
            # experiment_path / "deg_analysis" / "gene_stats_malignant_cibersortx.parquet"
        )
