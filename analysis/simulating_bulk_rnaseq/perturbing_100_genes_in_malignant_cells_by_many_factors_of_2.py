import logging

import cloudpathlib
import numpy as np
import pandas as pd

import helpers
from helpers import datasets
from helpers.simulating_bulk_rnaseq import creating_mixtures
from helpers.simulating_bulk_rnaseq.gene_perturbation import (
    determine_genes_to_perturb,
    perturb_scrnaseq_gene_expression_by_scaling_factor_in_cell_type,
)

logger = logging.getLogger(__name__)


def configure_logging():
    helpers.logging.configure_logging()
    logger.setLevel("DEBUG")
    logging.getLogger("helpers").setLevel("DEBUG")
    logging.getLogger("helpers.creating_mixtures").setLevel("INFO")


def load_scrnaseq_and_filter_genes():
    df_scrnaseq, df_sc_metadata = datasets.load_jerby_arnon(
        ref_genome="hg19", units="tpm"
    )
    df_bulkrnaseq_tcga_skcm = datasets.load_tcga_skcm_hg19_scaled_estimate_firebrowse()
    logger.debug("determining good genes")
    good_genes = helpers.data_io_and_formatting.qa_gene_filtering.get_good_genes(
        df_bulkrnaseq_tcga_skcm, df_scrnaseq, 0.5
    )
    logger.debug("limiting scRNA-seq data to good genes")
    df_scrnaseq = df_scrnaseq.loc[list(sorted(good_genes))]
    logger.debug("Shape of scRNA-seq data after filtering: %s", df_scrnaseq.shape)
    return df_scrnaseq, df_sc_metadata


if __name__ == "__main__":
    configure_logging()
    timestamp_str = helpers.useful_small_things.make_a_nice_timestamp_of_now()
    N = 50
    rng = np.random.default_rng(seed=0)

    # load data
    tcga_skcm_mets_fractions = (
        datasets.tcga_skcm.load_tcga_skcm_mets_fractions_from_csx()
    )
    df_sample_fractions = tcga_skcm_mets_fractions.sample(
        N, replace=False, random_state=rng
    )
    logger.debug("Loaded sample fractions, such as: %s", df_sample_fractions.head())
    df_scrnaseq, df_sc_metadata = load_scrnaseq_and_filter_genes()
    genes_to_perturb = determine_genes_to_perturb(df_scrnaseq, df_sc_metadata, rng)

    for scaling_factor in [0.125, 0.25, 0.5, 2.0, 4.0, 8.0]:
        df_scrnaseq_perturbed = (
            perturb_scrnaseq_gene_expression_by_scaling_factor_in_cell_type(
                df_scrnaseq,
                df_sc_metadata,
                "Malignant",
                genes_to_perturb,
                scaling_factor,
            )
        )

        # simulate bulk RNA-seq
        rng = np.random.default_rng(seed=0)
        (
            df_simulated_bulkrnaseq,
            simulated_cell_type_geps,
        ) = creating_mixtures.make_mixtures(
            df_scrnaseq_perturbed, df_sc_metadata, df_sample_fractions, rng=rng
        )

        # save stuff
        path_root = (
            cloudpathlib.CloudPath("gs://liulab/data/simulated/")
            / "perturbing_100_genes_in_malignant_cells_by_many_factors_of_2"
            / timestamp_str
            / f"scaling_factor={scaling_factor:.3f}"
        )
        df_simulated_bulkrnaseq.to_parquet(str(path_root / "bulk_rnaseq.parquet"))
        df_simulated_cell_type_geps = pd.concat(
            simulated_cell_type_geps, names=["sample_id"]
        )
        df_simulated_cell_type_geps.to_parquet(
            str(path_root / "cell_type_geps.parquet")
        )
        df_sample_fractions.to_parquet(str(path_root / "fractions.parquet"))
        perturbed_genes = pd.Series(genes_to_perturb, name="gene")
