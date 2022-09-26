import logging

import upath

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
    logging.getLogger("helpers.creating_mixtures").setLevel("INFO")

    # load data
    cohorts = {
        "control_1": upath.UPath(
            "gs://liulab/simulated/control/20220926_22h06m11s/"
        ),
        "control_2": upath.UPath(
            "gs://liulab/simulated/control/20220926_22h07m38s/"
        ),
    }
    df_bulk_rnaseq = load_and_concatenate_bulk_rnaseq(cohorts)
    assert df_bulk_rnaseq.shape[1] == 100, df_bulk_rnaseq.shape
    logger.debug(df_bulk_rnaseq)
    logger.debug(df_bulk_rnaseq.sum())
    df_fractions = load_and_concatenate_fractions(cohorts)
    assert df_fractions.shape[0] == 100, df_fractions.shape
    logger.debug(df_fractions)

    # define output path
    timestamp_str = helpers.useful_small_things.make_a_nice_timestamp_of_now()
    path_to_save_results_in_cloud = (
        upath.UPath("gs://liulab/cibersortx/negative_control") / timestamp_str
    )

    # run CIBERSORTx
    helpers.running_cibersortx.hires_only.run_and_upload_from_dataframes(
        str(path_to_save_results_in_cloud), df_bulk_rnaseq, df_fractions
    )
