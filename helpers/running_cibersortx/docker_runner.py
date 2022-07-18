import logging
import os
import pathlib
import tempfile

import dask.dataframe as dd
import docker
import pandas as pd
from cloudpathlib import AnyPath
from google.cloud import storage

import helpers

logger = logging.getLogger(__name__)


def create_csx_mixture_tsv(df_bulk_rnaseq: pd.DataFrame, uri_file: str) -> str:
    # tmp_dir = tmpfile.mkdtemp()
    # uri_csv_file = tmp_dir / "bulk_rnaseq.csv"
    # uri_csv_file = tempfile.mkstemp()
    df_bulk_rnaseq = df_bulk_rnaseq.rename_axis(index="Gene", columns=None)
    df_bulk_rnaseq *= 1e6 / df_bulk_rnaseq.sum()
    df_bulk_rnaseq = df_bulk_rnaseq.reset_index()
    df_bulk_rnaseq.to_csv(uri_file, index=False, sep="\t")
    return uri_file


def create_csx_refsample_tsv(
    df_refsample_sc_rnaseq: pd.DataFrame,
    df_refsample_sc_metadata: pd.DataFrame,
    uri_file: str,
) -> str:
    # tmp_dir = tempfile.mkdtemp()
    # uri_csv_file = os.path.join(tmp_dir, "refsample_sc_rnaseq.tsv")
    # uri_csv_file = tempfile.mkstemp()

    def look_up_cell_type(single_cell_id):
        return df_refsample_sc_metadata.loc[single_cell_id][helpers.columns.CELL_TYPE]

    df_refsample_sc_rnaseq = df_refsample_sc_rnaseq.rename(columns=look_up_cell_type)
    df_refsample_sc_rnaseq = df_refsample_sc_rnaseq.rename_axis(
        index="Gene", columns=None
    )
    df_refsample_sc_rnaseq = df_refsample_sc_rnaseq.reset_index()
    df_refsample_sc_rnaseq.to_csv(uri_file, index=False, sep="\t")
    return uri_file


def run_fractions_in_prepared_local_directory(csx_dir):
    run_kwargs = dict(
        user=os.getuid(),
        volumes=[
            f"{csx_dir}/in:/src/data",
            f"{csx_dir}/out:/src/outdir",
        ],
        auto_remove=True,
    )
    command_arguments = " \\ \n".join(
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
    message = client.containers.run(
        "cibersortx/fractions:latest", command_arguments, **run_kwargs
    )
    print(message.decode("utf-8"))
    # for message in container.logs(follow=True):
    #     print(message)
    # container.wait()


def run_fractions_and_upload(
    uri_bulk_rnaseq: str,
    uri_refsample_sc_rnaseq: str,
    uri_save_job_files_to: str,
):
    with tempfile.TemporaryDirectory() as tmp_dir:
        path = pathlib.Path(tmp_dir)
        (path / "in").mkdir()
        (path / "out").mkdir()
        AnyPath(uri_bulk_rnaseq).copy(path / "in" / "inputbulkrnaseq.tsv")
        AnyPath(uri_refsample_sc_rnaseq).copy(path / "in" / "inputrefscrnaseq.tsv")
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
