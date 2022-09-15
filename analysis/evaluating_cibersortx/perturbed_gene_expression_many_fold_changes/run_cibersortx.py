import logging

import cloudpathlib

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
    timestamp_str = helpers.useful_small_things.make_a_nice_timestamp_of_now()

    # load data
    path_root = cloudpathlib.CloudPath("gs://liulab/data/simulated/")
    path_unperturbed = path_root / "50_samples_no_perturbations" / "YYYY-MM-DD_HH:MM:SS"
    path_perturbed = (
        path_root
        / "perturbing_100_genes_in_malignant_cells_by_many_factors_of_2"
        / "YYYY-MM-DD_HH:MM:SS"
    )
    for scaling_factor in [0.125, 0.25, 0.5, 2.0, 4.0, 8.0]:
        simulated_cohorts = {
            "unperturbed": path_unperturbed,
            "perturbed": path_perturbed / f"scaling_factor={scaling_factor:.3f}",
        }
        df_bulk_rnaseq = load_and_concatenate_bulk_rnaseq(simulated_cohorts)
        assert df_bulk_rnaseq.shape[1] == 100, df_bulk_rnaseq.shape
        logger.debug(df_bulk_rnaseq)
        logger.debug(df_bulk_rnaseq.sum())
        df_fractions = load_and_concatenate_fractions(simulated_cohorts)
        assert df_fractions.shape[0] == 100, df_fractions.shape
        logger.debug(df_fractions)

        # define output path
        path_to_save_results_in_cloud = (
            cloudpathlib.CloudPath(
                "gs://liulab/evaluating_cibersortx/perturbing_100_genes_in_malignant_cells_by_many_factors_of_2"
            )
            / timestamp_str
            / f"scaling_factor={scaling_factor:.3f}",
        )

        # run CIBERSORTx
        helpers.running_cibersortx.hires_only.run_and_upload_from_dataframes(
            str(path_to_save_results_in_cloud), df_bulk_rnaseq, df_fractions
        )
