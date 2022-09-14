import logging
import tempfile

import cloudpathlib
import pandas as pd

import helpers

logger = logging.getLogger(__name__)


def load_and_concatenate_bulk_rnaseq(
    cohort_paths: dict[str, cloudpathlib.CloudPath]
) -> pd.DataFrame:
    df_simulated_bulkrnaseq = pd.concat(
        {
            cohort: pd.read_parquet(path / "bulk_rnaseq.parquet")
            for cohort, path in cohort_paths.items()
        },
        axis="columns",
    )
    df_simulated_bulkrnaseq.columns = (
        df_simulated_bulkrnaseq.columns.to_flat_index().map(lambda x: "/".join(x))
    )
    df_simulated_bulkrnaseq.rename_axis(index="gene_symbol", columns="sample_id", inplace=True)
    logger.debug("Shape of simulated data: %s", df_simulated_bulkrnaseq.shape)
    assert df_simulated_bulkrnaseq.shape[1] == 100, df_simulated_bulkrnaseq.shape
    return df_simulated_bulkrnaseq


def load_and_concatenate_fractions(
    cohort_paths: dict[str, cloudpathlib.CloudPath]
) -> pd.DataFrame:
    df_simulated_fractions = pd.concat(
        {
            cohort: pd.read_parquet(path / "fractions.parquet")
            for cohort, path in cohort_paths.items()
        },
        axis="rows",
    )
    df_simulated_fractions.index = df_simulated_fractions.index.to_flat_index().map(
        lambda x: "/".join(x)
    )
    df_simulated_fractions.rename_axis(index="sample_id", columns="cell_type", inplace=True)
    logger.debug("Shape of simulated data: %s", df_simulated_fractions.shape)
    assert df_simulated_fractions.shape[0] == 100, df_simulated_fractions.shape
    return df_simulated_fractions


if __name__ == "__main__":
    helpers.logging.configure_logging()
    logger.setLevel("DEBUG")
    logging.getLogger("helpers").setLevel("DEBUG")
    logging.getLogger("helpers.creating_mixtures").setLevel("INFO")

    # load simulated data
    simulated_cohorts = {
        "control": cloudpathlib.CloudPath(
            "gs://liulab/data/simulated/50_samples_no_perturbations/2022-09-13_21:37:53"
        ),
        "perturbed_2x": cloudpathlib.CloudPath(
            "gs://liulab/data/simulated/50_samples_100_genes_perturbed_2x_in_malignant_cells/2022-09-13_21:36:32"
        ),
    }
    df_simulated_bulk_rnaseq = load_and_concatenate_bulk_rnaseq(simulated_cohorts)
    logger.debug(df_simulated_bulk_rnaseq)
    df_fractions = load_and_concatenate_fractions(simulated_cohorts)
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
        str(path_to_save_results_in_cloud), df_simulated_bulk_rnaseq, df_fractions
    )
