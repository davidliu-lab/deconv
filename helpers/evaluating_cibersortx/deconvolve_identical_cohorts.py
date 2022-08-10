import cloudpathlib
import pandas as pd

from helpers import creating_mixtures, datasets, running_cibersortx

sc_data, sc_metadata = datasets.load_jerby_arnon(ref_genome="hg19", units="tpm")
tcga_skcm_fractions = datasets.load_tcga_skcm_fractions_from_csx()
tcga_skcm_metastatic_sample_barcodes = (
    datasets.load_tcga_skcm_metastatic_sample_barcodes()
)
tcga_skcm_fractions_mets = datasets.get_tcga_skcm_metastatic_sample_metadata()

"""
todo
- limit tcga skcm samples to metastases (where does filtering happen elsewhere?)
"""


def run_trial(
    sc_data: pd.DataFrame,
    sc_metadata: pd.DataFrame,
    fractions: pd.DataFrame,
    trial_id: int,
) -> None:
    sc_data = sc_data.sample(1000)
    path_results_prefix = cloudpathlib.AnyPath(
        f"gs://liulab/evaluating_cibersortx/identical_cohorts/{trial_id}"
    )
    df_pseudobulks = creating_mixtures.make_mixtures(
        sc_data, sc_metadata, tcga_skcm_fractions_mets
    )
    uri_pseudobulks = str(
        path_results_prefix / "csx_fractions" / "in" / "pseudobulks.tsv"
    )
    uri_refsample_sc_rnaseq = str(
        path_results_prefix / "csx_fractions" / "in" / "refscrnaseq.tsv"
    )
    uri_sigmatrix = str(
        path_results_prefix
        / "csx_fractions"
        / "outdir"
        / "CIBERSORTx_inputrefscrnaseq_inferred_phenoclasses.CIBERSORTx_refscrnaseq_inferred_refsample.bm.K999.txt"
    )
    running_cibersortx.creating_input_files.create_csx_mixtures_tsv(
        df_pseudobulks, uri_pseudobulks
    )
    # uri_results_fractions =
    # uri_results_hires = f"gs://liulab/evaluating_cibersortx/identical_cohorts/{trial_id}/results_hires"
    uri_results_hires = str(path_results_prefix / "csx_results" / "hires")
    running_cibersortx.sigmatrix_and_fractions.run_fractions_and_upload(
        uri_bulk_rnaseq,
        uri_refsample_sc_rnaseq,
        uri_save_job_files_to,
    )
    # uri_pseudobulks, uri_refsample_sc_rnaseq, uri_results_fractions
    running_cibersortx.hires_with_fractions.run_and_upload(
        uri_save_job_files_to,
        uri_bulk_rnaseq,
        uri_sigmatrix,
        uri_sourcegeps,
    )
