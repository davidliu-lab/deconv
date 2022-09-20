import logging

import pandas as pd
import upath

import helpers
from helpers import datasets
from helpers.simulating_bulk_rnaseq import creating_mixtures
from helpers.simulating_bulk_rnaseq.saving_to_cloud import save_simulated_data

logger = logging.getLogger(__name__)


def make_pseudobulks_and_run_hires(
    df_scrnaseq: pd.DataFrame,
    df_scrnaseq_metadata: pd.DataFrame,
    pseudobulk_sample_fractions: pd.DataFrame,
    path_to_save_results_in_cloud: upath.UPath,
) -> None:
    df_simulated_bulkrnaseq, simulated_cell_type_geps = creating_mixtures.make_mixtures(
        df_scrnaseq, df_scrnaseq_metadata, pseudobulk_sample_fractions
    )
    save_simulated_data(
        df_simulated_bulkrnaseq,
        simulated_cell_type_geps,
        path_to_save_results_in_cloud / "all_simulated_data",
    )
    path_bulk_rnaseq = path_to_save_results_in_cloud / "pseudobulk_rnaseq.tsv"
    helpers.running_cibersortx.creating_input_files.create_csx_mixtures_tsv(
        df_simulated_bulkrnaseq, path_bulk_rnaseq
    )
    helpers.running_cibersortx.hires_with_fractions.run_and_upload(
        str(path_to_save_results_in_cloud),
        path_bulk_rnaseq,
        uri_sigmatrix="gs://liulab/data/pseudobulk_evaluation/csx_runs/tcga_skcm/out/CIBERSORTx_inputrefscrnaseq_inferred_phenoclasses.CIBERSORTx_inputrefscrnaseq_inferred_refsample.bm.K999.txt",
        uri_sourcegeps="gs://liulab/data/pseudobulk_evaluation/csx_runs/tcga_skcm/out/CIBERSORTx_cell_type_sourceGEP.txt",
    )


if __name__ == "__main__":
    helpers.logging.configure_logging()
    logger.setLevel("DEBUG")
    logging.getLogger("helpers").setLevel("DEBUG")
    logging.getLogger("helpers.creating_mixtures").setLevel("INFO")

    logger.debug("loading pseudobulk sample metadata")
    tcga_skcm_fractions_mets = datasets.get_tcga_skcm_metastatic_sample_metadata()
    tcga_skcm_fractions = datasets.load_tcga_skcm_fractions_from_csx()
    pseudobulk_sample_fractions = tcga_skcm_fractions.loc[
        tcga_skcm_fractions_mets["aliquot_barcode"]
    ]
    logger.debug("number of pseudobulk samples: %s", len(pseudobulk_sample_fractions))
    df_scrnaseq, df_sc_metadata = datasets.load_jerby_arnon(
        ref_genome="hg19", units="tpm"
    )
    logger.debug("loading TCGA-SKCM bulk RNA-seq data")
    df_bulkrnaseq_tcga_skcm = datasets.load_tcga_skcm_hg19_scaled_estimate_firebrowse()
    logger.debug("determining high-quality genes")
    good_genes = helpers.data_io_and_formatting.qa_gene_filtering.get_good_genes(
        df_bulkrnaseq_tcga_skcm, df_scrnaseq, 0.5
    )
    logger.debug("limiting scRNA-seq data to high-quality genes")
    df_scrnaseq = df_scrnaseq.loc[list(sorted(good_genes))]
    logger.debug("Shape of scRNA-seq data: %s", df_scrnaseq.shape)
    timestamp_str = helpers.useful_small_things.make_a_nice_timestamp_of_now()
    for ith_trial in range(2):
        path_results = (
            upath.UPath("gs://liulab/evaluating_cibersortx/identical_cohorts")
            / timestamp_str
            / str(ith_trial)
        )
        make_pseudobulks_and_run_hires(
            df_scrnaseq, df_sc_metadata, pseudobulk_sample_fractions, path_results
        )
