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


def run_fractions_in_prepared_local_directory(csx_dir):
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
            "--single_cell TRUE",
            "--refsample inputrefscrnaseq.txt",
            "--mixture inputbulkrnaseq.txt",
            "--rmbatchBmode TRUE",
            "--verbose TRUE",
            "--absolute FALSE",
        ]
    )
    client = docker.from_env()
    container = client.containers.run(
        "cibersortx/fractions:latest", command_arguments, **run_kwargs
    )
    for message in container.logs(follow=True, stream=True):
        print(message.decode("utf-8"), end="")
    container.wait()


def set_up_fractions_dir(uri_bulk_rnaseq, uri_refsample_sc_rnaseq, tmp_dir):
    csx_dir = pathlib.Path(tmp_dir)
    (csx_dir / "data").mkdir()
    (csx_dir / "outdir").mkdir()
    copy_file_maybe_in_the_cloud_to_local_path(
        uri_bulk_rnaseq, csx_dir / "data" / "inputbulkrnaseq.txt"
    )
    copy_file_maybe_in_the_cloud_to_local_path(
        uri_refsample_sc_rnaseq, csx_dir / "data" / "inputrefscrnaseq.txt"
    )
    return csx_dir


def run_fractions_and_upload(
    uri_bulk_rnaseq: str,
    uri_refsample_sc_rnaseq: str,
    uri_save_job_files_to: str,
):
    with tempfile.TemporaryDirectory() as tmp_dir:
        logger.debug(f"running fractions in {tmp_dir}")
        set_up_fractions_dir(uri_bulk_rnaseq, uri_refsample_sc_rnaseq, tmp_dir)
        run_fractions_in_prepared_local_directory(tmp_dir)
        storage_client = storage.Client()
        bucket = storage_client.bucket("liulab")
        copy_local_directory_to_gcs(tmp_dir, bucket, uri_save_job_files_to)


if __name__ == "__main__":
    helpers.logging.configure_logging()
    logging.getLogger("gcsfs").setLevel("INFO")
    logging.getLogger("google.cloud.bigquery").setLevel("DEBUG")
    logging.getLogger("helpers").setLevel("DEBUG")
    logging.getLogger("pandas").setLevel("DEBUG")
    logging.getLogger("pyarrow").setLevel("DEBUG")
    logger.setLevel("DEBUG")
    logger.debug("test debug-level message")

    uri_refsample = "gs://liulab/data/pseudobulk_evaluation/csx_input_files/refsample_jerby_arnon_hg19_tpm.tsv"

    logger.debug("run cibersortx on tcga skcm bulk rna-seq")
    run_fractions_and_upload(
        uri_bulk_rnaseq="gs://liulab/data/pseudobulk_evaluation/csx_input_files/bulk_rnaseq_tcga_skcm.tsv",
        uri_refsample_sc_rnaseq=uri_refsample,
        uri_save_job_files_to="gs://liulab/data/pseudobulk_evaluation/csx_runs/tcga_skcm/",
    )

    logger.debug("run cibersortx on pseudobulk bulk rna-seq")
    run_fractions_and_upload(
        uri_bulk_rnaseq="gs://liulab/data/pseudobulk_evaluation/csx_input_files/bulk_rnaseq_pseudobulk.tsv",
        uri_refsample_sc_rnaseq=uri_refsample,
        uri_save_job_files_to="gs://liulab/data/pseudobulk_evaluation/csx_runs/pseudobulks/n_cells=5/malignant_from_one_sample=True/",
    )
