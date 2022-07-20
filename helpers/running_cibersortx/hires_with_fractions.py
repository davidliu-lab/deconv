import tempfile

from .copying_to_gcs import copy_local_directory_to_gcs


def run_and_upload(
    uri_save_job_files_to, uri_bulk_rnaseq, uri_sigmatrix, uri_sourcegeps
):
    with tempfile.TemporaryDirectory() as tmp_dir:
        set_up_csx_dir(tmp_dir, uri_bulk_rnaseq, uri_fractions, uri_refsample)
        run(tmp_dir)
        copy_local_directory_to_gcs(tmp_dir, uri_save_job_files_to)


def set_up_csx_dir(csx_dir, uri_bulk_rnaseq, uri_sigmatrix, uri_sourcegeps):
    csx_path = pathlib.Path(csx_dir)
    (csx_path / "data").mkdir()
    (csx_path / "outdir").mkdir()
    copy_file_maybe_in_the_cloud_to_local_path(
        uri_bulk_rnaseq, csx_dir / "data" / "bulkrnaseq.txt"
    )
    copy_file_maybe_in_the_cloud_to_local_path(
        uri_sigmatrix, csx_dir / "data" / "sigmatrix.txt"
    )
    copy_file_maybe_in_the_cloud_to_local_path(
        uri_sourcegeps, csx_dir / "data" / "sourcegeps.txt"
    )


def run(csx_dir):
    run_kwargs = dict(
        auto_remove=True,
        detach=True,
        user=f"{os.getuid()}:{os.getgid()}",
        volumes=[
            f"{csx_dir}/data:/src/data",
            f"{csx_dir}/outdir:/src/outdir",
        ],
    )
    command_arguments = " ".join(
        [
            "--username lyronctk@stanford.edu",
            "--token dfeba2c8b9d61daebee5fa87026b8e56",
            "--mixture bulkrnaseq.txt",  # <file_name>  Mixture matrix [required]
            "--rmbatchBmode TRUE",  # only relevant if hires is also running fractions
            # "--verbose TRUE",  # not sure this is even an option?
            "--sourceGEPs sourcegeps.txt",  # <file_name>  Signature matrix GEPs for batch correction [default: sigmatrix]
            # "--groundtruth ",   # <file_name>  Ground truth GEPs [same labels as classes] [optional; default: none]            "--mixture bulkrnaseq.tsv",
            "--sigmatrix sigmatrix.txt",  # <file_name>  Signature matrix [required]
            # "--cibresults cibresults1.txt",
            # "--degclasses ",  # <file_name>  Run on two classes, specified by 1, 2, 0=skip [default: none]
            # "--filtered ", # <file_name>  Filtered GEPs from CIBERSORTxGEP [default: none]
            # "--window ",  # <int>   Window size for deconvolution [default: No. of cell types x 4]
            # "--useadjustedmixtures ",  # <bool>  If doing B-mode batch correction, use adjusted mixtures for GEP imputation [default: FALSE]
        ]
    )
    pass


if __name__ == "__main__":
    helpers.logging.configure_logging()
    logging.getLogger("gcsfs").setLevel("INFO")
    logging.getLogger("google.cloud.bigquery").setLevel("DEBUG")
    logging.getLogger("helpers").setLevel("DEBUG")
    logging.getLogger("pandas").setLevel("DEBUG")
    logging.getLogger("pyarrow").setLevel("DEBUG")
    logger.setLevel("DEBUG")

    logger.debug("run cibersortx hires (with fractions) on tcga skcm")
    run_upload(
        uri_save_job_files_to="gs://liulab/data/pseudobulk_evaluation/csx_runs/hires_with_fractions/tcga_skcm/",
        uri_bulk_rnaseq="gs://liulab/data/pseudobulk_evaluation/csx_input_files/bulk_rnaseq_tcga_skcm.tsv",
        uri_sigmatrix="gs://liulab/data/pseudobulk_evaluation/csx_runs/tcga_skcm/out/CIBERSORTx_inputrefscrnaseq_inferred_phenoclasses.CIBERSORTx_inputrefscrnaseq_inferred_refsample.bm.K999.txt",
        uri_sourcegeps="gs://liulab/data/pseudobulk_evaluation/csx_runs/tcga_skcm/out/CIBERSORTx_cell_type_sourceGEP.txt",
    )

    logger.debug("run cibersortx hires (with fractions) on pseudobulk bulk rna-seq")
    run_upload(
        uri_save_job_files_to="gs://liulab/data/pseudobulk_evaluation/csx_runs/hires_with_fractions/pseudobulks/n_cells=5/malignant_from_one_sample=True/",
        # uri_save_job_files_to="gs://liulab/data/pseudobulk_evaluation/csx_runs/hires_with_fractions/tcga_skcm/pseudobulks/n_cells=5/malignant_from_one_sample=True/",
        # uri_refsample_sc_rnaseq=uri_refsample,
        uri_bulk_rnaseq="gs://liulab/data/pseudobulk_evaluation/csx_input_files/bulk_rnaseq_pseudobulk.tsv",
        uri_sigmatrix="gs://liulab/data/pseudobulk_evaluation/csx_runs/pseudobulks/n_cells=5/malignant_from_one_sample=True/out/CIBERSORTx_inputrefscrnaseq_inferred_phenoclasses.CIBERSORTx_inputrefscrnaseq_inferred_refsample.bm.K999.txt",
        uri_sourcegeps="gs://liulab/data/pseudobulk_evaluation/csx_runs/pseudobulks/n_cells=5/malignant_from_one_sample=True/out/CIBERSORTx_cell_type_sourceGEP.txt",
    )
