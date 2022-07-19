import logging
import os
import pathlib
import shutil
import tempfile

import dask.dataframe as dd
import docker
import pandas as pd
from cloudpathlib import AnyPath
from google.cloud import storage

import helpers
from helpers.running_cibersortx.copying_to_gcs import (
    copy_local_directory_to_gcs,
    copy_file_maybe_in_the_cloud_to_local_path,
)
from helpers.running_cibersortx.creating_input_files import (
    create_csx_mixture_tsv,
    create_csx_refsample_tsv,
)

logger = logging.getLogger(__name__)


def run_fractions_in_prepared_local_directory(csx_dir):
    run_kwargs = dict(
        auto_remove=True,
        detach=True,
        user=f"{os.getuid()}:{os.getgid()}",
        volumes=[
            f"{csx_dir}/in:/src/data",
            f"{csx_dir}/out:/src/outdir",
        ],
    )
    command_arguments = " ".join(
        [
            "--username lyronctk@stanford.edu",
            "--token dfeba2c8b9d61daebee5fa87026b8e56",
            "--single_cell TRUE",
            "--refsample inputrefscrnaseq.tsv",
            "--mixture inputbulkrnaseq.tsv",
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


def run_fractions_and_upload(
    uri_bulk_rnaseq: str,
    uri_refsample_sc_rnaseq: str,
    uri_save_job_files_to: str,
):
    with tempfile.TemporaryDirectory() as tmp_dir:
        logger.debug(f"running cibersortx from {tmp_dir}")
        path = pathlib.Path(tmp_dir)
        (path / "in").mkdir()
        (path / "out").mkdir()
        copy_file_maybe_in_the_cloud_to_local_path(
            uri_bulk_rnaseq, path / "in" / "inputbulkrnaseq.tsv"
        )
        copy_file_maybe_in_the_cloud_to_local_path(
            uri_refsample_sc_rnaseq, path / "in" / "inputrefscrnaseq.tsv"
        )
        run_fractions_in_prepared_local_directory(tmp_dir)
        storage_client = storage.Client()
        bucket = storage_client.bucket("liulab")
        copy_local_directory_to_gcs(tmp_dir, bucket, uri_save_job_files_to)


if __name__ == "__main__":
    logging.getLogger("gcsfs").setLevel("INFO")
    logging.getLogger("google.cloud.bigquery").setLevel("DEBUG")
    logging.getLogger("helpers").setLevel("DEBUG")
    logging.getLogger("helpers.pseudobulk_evaluation.deg_analysis").setLevel("INFO")
    logging.getLogger("pandas").setLevel("DEBUG")
    logging.getLogger("pyarrow").setLevel("DEBUG")
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        "%(asctime)s %(process)d/%(threadName)s %(name)s %(levelname)s\n%(message)s"
    )
    handler.setFormatter(formatter)
    logging.getLogger().handlers = [handler]
    logger.setLevel("DEBUG")
    logger.debug("test debug-level message")
    # logger.debug("loading data (bulk and sc data)")
    # df_bulk_rnaseq_hg19_tpm = (
    #     helpers.datasets.load_tcga_skcm_hg19_scaled_estimate_firebrowse()
    # )
    # # df_bulk_rnaseq_hg19_tpm *= 1e6 / df_bulk_rnaseq_hg19_tpm.sum()
    # df_sc_rnaseq, df_sc_metadata = helpers.datasets.load_jerby_arnon_hg19_tpm()
    # # df_sc_rnaseq *= 1e6 / df_sc_rnaseq.sum()
    # logger.debug("creating bulk rna-seq tsv formatted for cibersortx")
    # create_csx_mixture_tsv(df_bulk_rnaseq_hg19_tpm, uri_csx_bulk_rnaseq)
    # logger.debug("creating refsample sc rna-seq tsv formatted for cibersortx")
    # create_csx_refsample_tsv(df_sc_rnaseq, df_sc_metadata, uri_csx_refsample)
    # uri_csx_bulk_rnaseq = "/tmp/csx_bulk_rnaseq.tsv"
    # uri_csx_refsample = "/tmp/csx_refsample.tsv"
    # logger.debug("running with docker")
    # run_high_resolution_and_upload(
    #     uri_csx_bulk_rnaseq,
    #     uri_csx_refsample,
    #     "gs://liulab/tmp/csx_job_files",
    # )
    run_fractions_in_prepared_local_directory("/home/jupyter/deconv/tmp")
