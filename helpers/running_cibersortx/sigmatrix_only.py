import logging
import pathlib
import tempfile

import docker
import gcsfs
import upath

import helpers

from .copying_to_gcs import copy_file_maybe_in_the_cloud_to_local_path

logger = logging.getLogger(__name__)


def set_up_local_dir(uri_refsample_sc_rnaseq: str, tmp_dir: str):
    path_csx = pathlib.Path(tmp_dir)
    (path_csx / "data").mkdir()
    (path_csx / "outdir").mkdir()
    copy_file_maybe_in_the_cloud_to_local_path(
        uri_refsample_sc_rnaseq, path_csx / "data" / "inputrefscrnaseq.txt"
    )


def run(csx_dir: str) -> None:
    """
    Run a cibersortx container in the given directory to generate a signature matrix

    :param csx_dir: a local directory with a "reference sample" file of scRNA-seq
    """
    logger.debug("generating signature matrix in %s", csx_dir)
    run_kwargs = dict(
        auto_remove=True,
        detach=True,
        # user=f"{os.getuid()}:{os.getgid()}",
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
            "--verbose TRUE",
        ]
    )
    client = docker.from_env()
    container = client.containers.run(
        "cibersortx/fractions:latest", command_arguments, **run_kwargs
    )
    for message in container.logs(follow=True, stream=True):
        print(message.decode("utf-8"), end="")
    container.wait()


def set_up_run_and_upload(
    uri_refsample_sc_rnaseq: str,
    uri_save_job_files_to: str,
):
    with tempfile.TemporaryDirectory() as tmp_dir:
        logger.debug("generating signature matrix in %s", tmp_dir)
        set_up_local_dir(uri_refsample_sc_rnaseq, tmp_dir)
        run(tmp_dir)
        gcsfs.GCSFileSystem().put(tmp_dir, uri_save_job_files_to, recursive=True)


if __name__ == "__main__":
    helpers.logging.configure_logging()
    logging.getLogger("helpers").setLevel("DEBUG")
    logger.setLevel("DEBUG")

    # uri_refsample = "gs://liulab/data/pseudobulk_evaluation/csx_input_files/refsample_jerby_arnon_hg19_tpm.tsv"
    base_path = (
        upath.UPath("gs://liulab/tmp") / helpers.useful_small_things.make_a_nice_timestamp_of_now()
    )
    uri_refsample = str(base_path / "sigmatrix_only_refsample.tsv")
    sc_data, sc_metadata = helpers.datasets.load_jerby_arnon()
    # sc_data = sc_data.iloc[::10, ::10]
    logger.debug("sc_data has shape %s", sc_data.shape)
    for i in range(10):
        logger.debug("making sigmat: trial %d", i)
        sc_data_i = sc_data.iloc[:, i::10]
        logger.debug("shape of sc_data_i: %s", sc_data_i.shape)
        helpers.running_cibersortx.creating_input_files.create_csx_refsample_tsv(
            sc_data_i, sc_metadata, uri_refsample
        )
        set_up_run_and_upload(uri_refsample, str(base_path / str(i) / "csx_outputs"))
