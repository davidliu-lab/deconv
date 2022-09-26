import logging

import numpy as np
import pandas as pd
import upath

import helpers
from helpers import datasets
from helpers.simulating_bulk_rnaseq import creating_mixtures

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    helpers.logging.configure_logging()
    logger.setLevel("DEBUG")
    logging.getLogger("helpers").setLevel("DEBUG")
    logging.getLogger("helpers.creating_mixtures").setLevel("INFO")
    timestamp_str = helpers.useful_small_things.make_a_nice_timestamp_of_now()
    N = 50
    rng = np.random.default_rng(seed=0)

    # load data
    df_fractions_tcga_skcm_mets = datasets.tcga_skcm.load_fractions_mets_only()
    df_scrnaseq, df_sc_metadata = datasets.jerby_arnon.load_scrnaseq_and_filter_genes()

    # randomly sample fractions
    df_fractions = creating_mixtures.randomly_sample_fractions(
        df_fractions_tcga_skcm_mets, N, rng
    )
    logger.debug("Randomly sampled fractions, such as: %s", df_fractions.head())

    # simulate bulk RNA-seq
    df_bulk_rnaseq, cell_type_geps = creating_mixtures.make_mixtures(
        df_scrnaseq, df_sc_metadata, df_fractions, rng=rng
    )
    df_cell_type_geps = pd.concat(cell_type_geps, names=["sample_id"])

    # save stuff
    path_root = upath.UPath("gs://liulab/simulated/control") / timestamp_str
    logger.debug("saving data to %s", path_root)
    logger.debug("saving df_bulk_rnaseq")
    df_bulk_rnaseq.to_parquet(path_root / "bulk_rnaseq.parquet")
    logger.debug("saving df_cell_type_geps")
    df_cell_type_geps.to_parquet(path_root / "cell_type_geps.parquet")
    logger.debug("saving df_fractions")
    df_fractions.to_parquet(path_root / "fractions.parquet")
