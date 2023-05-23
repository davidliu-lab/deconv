import logging

from upath import UPath

import helpers
from helpers.data_io_and_formatting.concatenating import (
    load_and_concatenate_bulk_rnaseq,
    load_and_concatenate_fractions,
)

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    helpers.logging.configure_logging()
    logger.setLevel("DEBUG")
    logging.getLogger("helpers").setLevel("DEBUG")
    timestamp_str = helpers.useful_small_things.make_a_nice_timestamp_of_now()

    # load data
    simulated_datasets = {
        "control_1": UPath("gs://liulab/simulated/control/20220927_15h06m39s/seed=0"),
        "control_2": UPath("gs://liulab/simulated/control/20220927_15h21m19s/seed=1"),
    }
    df_bulk_rnaseq = load_and_concatenate_bulk_rnaseq(simulated_datasets)
    assert df_bulk_rnaseq.shape[1] == 100, df_bulk_rnaseq.shape
    logger.debug(df_bulk_rnaseq)
    logger.debug(df_bulk_rnaseq.sum())
    df_fractions = load_and_concatenate_fractions(simulated_datasets)
    assert df_fractions.shape[0] == 100, df_fractions.shape
    logger.debug(df_fractions)

    # define output path
    path_results = UPath("gs://liulab/cibersortx/negative_control") / timestamp_str

    # run CIBERSORTx
    helpers.running_cibersortx.hires_only.run_and_upload_from_dataframes(
        df_bulk_rnaseq, df_fractions, path_results
    )
