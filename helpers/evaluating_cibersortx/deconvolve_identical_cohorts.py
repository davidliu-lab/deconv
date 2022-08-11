import datetime
import logging

import cloudpathlib
import pandas as pd

import helpers
from helpers import creating_mixtures, datasets, running_cibersortx

logger = logging.getLogger(__name__)


def make_pseudobulks_and_run_hires(
    sc_data: pd.DataFrame,
    sc_metadata: pd.DataFrame,
    pseudobulk_sample_fractions: pd.DataFrame,
    path_results: cloudpathlib.AnyPath,
) -> None:
    df_pseudobulk_rnaseq, cell_type_geps = creating_mixtures.make_mixtures(
        sc_data, sc_metadata, pseudobulk_sample_fractions
    )
    uri_bulk_rnaseq = str(path_results / "pseudobulk_rnaseq.tsv")
    helpers.running_cibersortx.creating_input_files.create_csx_mixtures_tsv(
        df_pseudobulk_rnaseq, uri_bulk_rnaseq
    )
    helpers.running_cibersortx.hires_with_fractions.run_and_upload(
        str(path_results),
        uri_bulk_rnaseq,
        uri_sigmatrix="gs://liulab/data/pseudobulk_evaluation/csx_runs/tcga_skcm/out/CIBERSORTx_inputrefscrnaseq_inferred_phenoclasses.CIBERSORTx_inputrefscrnaseq_inferred_refsample.bm.K999.txt",
        uri_sourcegeps="gs://liulab/data/pseudobulk_evaluation/csx_runs/tcga_skcm/out/CIBERSORTx_cell_type_sourceGEP.txt",
    )


if __name__ == "__main__":
    helpers.logging.configure_logging()
    logger.setLevel("DEBUG")
    logging.getLogger("helpers").setLevel("DEBUG")

    logger.debug("loading pseudobulk sample metadata")
    tcga_skcm_fractions_mets = datasets.get_tcga_skcm_metastatic_sample_metadata()
    tcga_skcm_fractions = datasets.load_tcga_skcm_fractions_from_csx()
    pseudobulk_sample_fractions = tcga_skcm_fractions.loc[
        tcga_skcm_fractions_mets["aliquot_barcode"]
    ]
    logger.debug(f"number of pseudobulk samples: {len(pseudobulk_sample_fractions)}")
    sc_data, sc_metadata = datasets.load_jerby_arnon(ref_genome="hg19", units="tpm")
    timestamp_str = helpers.useful_small_things.make_a_nice_timestamp_of_now()
    for ith_trial in range(2):
        path_results = (
            cloudpathlib.AnyPath("gs://liulab/evaluating_cibersortx/identical_cohorts")
            / timestamp_str
            / str(ith_trial)
        )
        make_pseudobulks_and_run_hires(
            sc_data, sc_metadata, pseudobulk_sample_fractions, path_results
        )
