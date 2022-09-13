import logging
import tempfile

import cloudpathlib
import dask.dataframe as dd
import pandas as pd

import helpers

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    helpers.logging.configure_logging()
    logger.setLevel("DEBUG")
    logging.getLogger("helpers").setLevel("DEBUG")
    logging.getLogger("helpers.creating_mixtures").setLevel("INFO")

    # load simulated data
    uris = {
        "control": "gs://liulab/data/simulated/50_samples_no_perturbations/0/bulk_rnaseq.parquet",
        "perturbed_2x": "gs://liulab/data/simulated/50_samples_100_genes_perturbed_in_malignant_cells/0/bulk_rnaseq.parquet",
    }
    # df_simulated_bulkrnaseq = pd.concat([pd.read_parquet(uri) for uri in uris])
    df_simulated_bulkrnaseq = pd.concat(
        {cohort: pd.read_parquet(uri) for cohort, uri in uris.items()}, axis="columns"
    )
    df_simulated_bulkrnaseq.columns = (
        df_simulated_bulkrnaseq.columns.to_flat_index().map(lambda x: "/".join(x))
    )
    logger.debug("Shape of simulated data: %s", df_simulated_bulkrnaseq.shape)
    assert len(df_simulated_bulkrnaseq) == 100

    # define output path
    timestamp_str = helpers.useful_small_things.make_a_nice_timestamp_of_now()
    path_to_save_results_in_cloud = (
        cloudpathlib.CloudPath(
            "gs://liulab/evaluating_cibersortx/perturbed_gene_expression"
        )
        / timestamp_str
    )

    with tempfile.NamedTemporaryFile(suffix=".tsv") as tmp_file:
        uri_bulk_rnaseq = tmp_file.name
        helpers.running_cibersortx.creating_input_files.create_csx_mixtures_tsv(
            df_simulated_bulkrnaseq, uri_bulk_rnaseq
        )
        helpers.running_cibersortx.hires_with_fractions.run_and_upload(
            str(path_to_save_results_in_cloud),
            uri_bulk_rnaseq,
            uri_sigmatrix="gs://liulab/data/pseudobulk_evaluation/csx_runs/tcga_skcm/out/CIBERSORTx_inputrefscrnaseq_inferred_phenoclasses.CIBERSORTx_inputrefscrnaseq_inferred_refsample.bm.K999.txt",
            uri_sourcegeps="gs://liulab/data/pseudobulk_evaluation/csx_runs/tcga_skcm/out/CIBERSORTx_cell_type_sourceGEP.txt",
        )
