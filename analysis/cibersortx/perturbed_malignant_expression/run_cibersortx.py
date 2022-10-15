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
    # logging.getLogger("gcsfs").setLevel("DEBUG")
    logging.getLogger("upath").setLevel("DEBUG")
    logging.getLogger("helpers").setLevel("DEBUG")
    timestamp_str = helpers.useful_small_things.make_a_nice_timestamp_of_now()

    # load data
    path_control = UPath("gs://liulab/simulated/control/20220927_15h06m39s/seed=0")
    root_perterbed = UPath(
        "gs://liulab/simulated/perturbed_malignant_expression/20221015_21h52m40s"
    )
    paths_perturbed = [p.parent for p in root_perterbed.glob("**/bulk_rnaseq.parquet")]
    logger.debug(paths_perturbed)
    for path_perturbed in paths_perturbed:
        logger.debug("starting %s", path_perturbed)
        scaling_factor_str = path_perturbed.name
        assert scaling_factor_str[:8] == "log2_fc=", scaling_factor_str
        simulated_datasets = {
            "control": path_control,
            scaling_factor_str: path_perturbed,
        }
        logger.debug(simulated_datasets)
        df_bulk_rnaseq = load_and_concatenate_bulk_rnaseq(simulated_datasets)
        assert df_bulk_rnaseq.shape[1] == 100, df_bulk_rnaseq.shape
        logger.debug(df_bulk_rnaseq)
        logger.debug(df_bulk_rnaseq.sum())
        df_fractions = load_and_concatenate_fractions(simulated_datasets)
        assert df_fractions.shape[0] == 100, df_fractions.shape
        logger.debug(df_fractions)

        # define output path
        path_results = (
            UPath("gs://liulab/cibersortx/perturbed_malignant_expression")
            / timestamp_str
            / scaling_factor_str
        )
        logger.debug(path_results)

        # run CIBERSORTx
        helpers.running_cibersortx.hires_only.run_and_upload_from_dataframes(
            df_bulk_rnaseq, df_fractions, path_results
        )
