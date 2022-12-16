import logging

import numpy as np
import pandas as pd
from upath import UPath

import helpers
from helpers import datasets
from helpers.simulating_bulk_rnaseq import creating_mixtures
from helpers.simulating_bulk_rnaseq.gene_perturbation import (
    select_100_genes_at_least_somewhat_expressed_in_malignant,
    perturb_scrnaseq_gene_expression,
)

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    helpers.logging.configure_logging()
    logger.setLevel("DEBUG")
    logging.getLogger("helpers").setLevel("DEBUG")
    logging.getLogger("gcsfs").setLevel("INFO")
    logging.getLogger("helpers.creating_mixtures").setLevel("INFO")
    timestamp_str = helpers.useful_small_things.make_a_nice_timestamp_of_now()
    N = 50
    path_root = UPath("gs://liulab/simulated/perturbed_malignant_expression") / timestamp_str
    # load data
    df_fractions_tcga_skcm_mets = datasets.tcga_skcm.load_fractions_mets_only()
    df_scrnaseq, df_sc_metadata = datasets.jerby_arnon.load_scrnaseq_and_filter_genes()

    # randomly sample genes to perturb
    rng = np.random.default_rng(seed=0)
    genes_to_perturb = select_100_genes_at_least_somewhat_expressed_in_malignant(
        df_scrnaseq, df_sc_metadata, rng
    )

    # save genes to perturb
    logger.debug("saving genes_to_perturb")
    pd.DataFrame(genes_to_perturb).to_csv(path_root / "genes_perturbed.csv", index=False)

    for i, log2_fc in enumerate(
        [-3, -2, -1.5, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3]
    ):
        scaling_factor = 2**log2_fc
        seed = 10 + i
        rng = np.random.default_rng(seed=seed)
        path = path_root / f"seed={seed}" / f"log2_fc={log2_fc:.3f}"
        logger.debug("saving data to %s", path)

        # perturb gene expression
        df_scrnaseq_perturbed = perturb_scrnaseq_gene_expression(
            df_scrnaseq, df_sc_metadata, "Malignant", genes_to_perturb, scaling_factor
        )

        # randomly sample fractions
        df_fractions = df_fractions_tcga_skcm_mets.sample(N, replace=True, random_state=rng)
        sample_index = pd.Index([f"sample_{i:03d}" for i in range(N)], name="sample_id")
        df_fractions.set_index(sample_index, inplace=True)
        logger.debug("Randomly sampled fractions, such as: %s", df_fractions.head())
        logger.debug("saving df_fractions")
        df_fractions.to_parquet(path / "fractions.parquet")

        # simulate bulk RNA-seq
        df_bulk_rnaseq, cell_type_geps = creating_mixtures.make_mixtures(
            df_scrnaseq_perturbed, df_sc_metadata, df_fractions, rng=rng
        )
        df_cell_type_geps = pd.concat(cell_type_geps, names=["sample_id"])
        logger.debug("saving df_bulk_rnaseq")
        df_bulk_rnaseq.to_parquet(path / "bulk_rnaseq.parquet")
        logger.debug("saving df_cell_type_geps")
        df_cell_type_geps.to_parquet(path / "cell_type_geps.parquet")
