import json
import logging
import os
import pathlib
import tempfile

import docker
from google.cloud import storage

import helpers
from helpers.running_cibersortx.copying_to_gcs import (
    copy_file_maybe_in_the_cloud_to_local_path,
    copy_local_directory_to_gcs,
)

logger = logging.getLogger(__name__)


def run_and_upload(
    uri_save_job_files_to, uri_bulk_rnaseq, uri_sigmatrix, uri_sourcegeps
):
    with tempfile.TemporaryDirectory(
        # dir="/Users/william/src/deconv/tmp/",
    ) as tmp_dir:
        logger.debug(f"tmp_dir: {tmp_dir}")
        logger.debug(f"watch -n 0.1 tree -ghpu {tmp_dir}")
        set_up_csx_dir(tmp_dir, uri_bulk_rnaseq, uri_sigmatrix, uri_sourcegeps)
        run(tmp_dir)
        storage_client = storage.Client()
        bucket = storage_client.bucket("liulab")
        copy_local_directory_to_gcs(tmp_dir, bucket, uri_save_job_files_to)


def set_up_csx_dir(csx_dir, uri_bulk_rnaseq, uri_sigmatrix, uri_sourcegeps):
    csx_path = pathlib.Path(csx_dir)
    (csx_path / "data").mkdir()
    (csx_path / "outdir").mkdir()
    copy_file_maybe_in_the_cloud_to_local_path(
        uri_bulk_rnaseq, csx_path / "data" / "bulkrnaseq.txt"
    )
    copy_file_maybe_in_the_cloud_to_local_path(
        uri_sigmatrix, csx_path / "data" / "sigmatrix.txt"
    )
    copy_file_maybe_in_the_cloud_to_local_path(
        uri_sourcegeps, csx_path / "data" / "sourcegeps.txt"
    )


def run(csx_dir):
    run_kwargs = dict(
        auto_remove=True,
        detach=True,
        # user=f"{os.getuid()}:{os.getgid()}",
        volumes=[
            f"{csx_dir}/data:/src/data",
            f"{csx_dir}/outdir:/src/outdir",
        ],
    )
    logger.debug(f"run_kwargs:\n{json.dumps(run_kwargs, indent=2, sort_keys=True)}")
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
    client = docker.from_env()
    container = client.containers.run(
        "cibersortx/hires:latest", command_arguments, **run_kwargs
    )
    for message in container.logs(follow=True, stream=True):
        print(message.decode("utf-8"), end="")
    container.wait()


if __name__ == "__main__":
    helpers.logging.configure_logging()
    logging.getLogger("gcsfs").setLevel("INFO")
    logging.getLogger("google.cloud.bigquery").setLevel("DEBUG")
    logging.getLogger("helpers").setLevel("DEBUG")
    logging.getLogger("pandas").setLevel("DEBUG")
    logging.getLogger("pyarrow").setLevel("DEBUG")
    logger.setLevel("DEBUG")

    data_path_prefix = "gs://liulab/data"
    # data_path_prefix = "/Users/william/Downloads/liulab_mirror/data"

    logger.debug("run cibersortx hires (with fractions) on tcga skcm")
    run_and_upload(
        uri_save_job_files_to=f"{data_path_prefix}/pseudobulk_evaluation/csx_runs/hires_with_fractions/tcga_skcm/",
        uri_bulk_rnaseq=f"{data_path_prefix}/pseudobulk_evaluation/csx_input_files/bulk_rnaseq_tcga_skcm.tsv",
        uri_sigmatrix=f"{data_path_prefix}/pseudobulk_evaluation/csx_runs/tcga_skcm/out/CIBERSORTx_inputrefscrnaseq_inferred_phenoclasses.CIBERSORTx_inputrefscrnaseq_inferred_refsample.bm.K999.txt",
        uri_sourcegeps=f"{data_path_prefix}/pseudobulk_evaluation/csx_runs/tcga_skcm/out/CIBERSORTx_cell_type_sourceGEP.txt",
    )

    logger.debug("run cibersortx hires (with fractions) on pseudobulk bulk rna-seq")
    run_and_upload(
        uri_save_job_files_to=f"{data_path_prefix}/pseudobulk_evaluation/csx_runs/hires_with_fractions/pseudobulks/n_cells=5/malignant_from_one_sample=True/",
        uri_bulk_rnaseq=f"{data_path_prefix}/pseudobulk_evaluation/csx_input_files/bulk_rnaseq_pseudobulk.tsv",
        uri_sigmatrix=f"{data_path_prefix}/pseudobulk_evaluation/csx_runs/pseudobulks/n_cells=5/malignant_from_one_sample=True/out/CIBERSORTx_inputrefscrnaseq_inferred_phenoclasses.CIBERSORTx_inputrefscrnaseq_inferred_refsample.bm.K999.txt",
        uri_sourcegeps=f"{data_path_prefix}/pseudobulk_evaluation/csx_runs/pseudobulks/n_cells=5/malignant_from_one_sample=True/out/CIBERSORTx_cell_type_sourceGEP.txt",
    )
