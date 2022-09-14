import logging

import cloudpathlib

import helpers
from helpers.data_formatting.concatenating import (
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
    simulated_cohorts = {
        "control": cloudpathlib.CloudPath(
            "gs://liulab/data/simulated/50_samples_no_perturbations/2022-09-13_21:37:53"
        ),
        "perturbed_2x": cloudpathlib.CloudPath(
            "gs://liulab/data/simulated/50_samples_100_genes_perturbed_2x_in_malignant_cells/2022-09-13_21:36:32"
        ),
    }
    df_bulk_rnaseq = load_and_concatenate_bulk_rnaseq(simulated_cohorts)
    assert df_bulk_rnaseq.shape[1] == 100, df_bulk_rnaseq.shape
    logger.debug(df_bulk_rnaseq)
    logger.debug(df_bulk_rnaseq.sum())
    df_fractions = load_and_concatenate_fractions(simulated_cohorts)
    assert df_fractions.shape[0] == 100, df_fractions.shape
    logger.debug(df_fractions)

    # define output path
    timestamp_str = helpers.useful_small_things.make_a_nice_timestamp_of_now()
    path_to_save_results_in_cloud = (
        cloudpathlib.CloudPath(
            "gs://liulab/evaluating_cibersortx/perturbed_gene_expression/2x"
        )
        / timestamp_str
    )

    # run CIBERSORTx
    helpers.running_cibersortx.hires_only.run_and_upload_from_dataframes(
        str(path_to_save_results_in_cloud), df_bulk_rnaseq, df_fractions
    )
