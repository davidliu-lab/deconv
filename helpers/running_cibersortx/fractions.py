import logging
import os
import pathlib
import tempfile

import docker
import pandas as pd
from upath import UPath

from helpers.running_cibersortx.creating_input_files import create_csx_mixtures_tsv

logger = logging.getLogger(__name__)


def run_and_upload_from_dataframe(
    bulk_rnaseq: pd.DataFrame,
    sigmatrix_path: UPath,
    sourcegeps_path: UPath,
    path_target_on_gcs: UPath,
):
    with tempfile.TemporaryDirectory() as tmp_dir:
        logger.debug("tmp_dir: %s", tmp_dir)
        csx_path = pathlib.Path(tmp_dir)
        (csx_path / "data").mkdir()
        (csx_path / "outdir").mkdir()
        create_csx_mixtures_tsv(bulk_rnaseq, csx_path / "data" / "bulkrnaseq.txt")
        sigmatrix_path.fs.get(csx_path / "outdir" / "sigmatrix.txt")
        sourcegeps_path.fs.get(csx_path / "outdir" / "sourcegeps.txt")
        run(tmp_dir)
        logger.debug("copying tmp_dir %s to %s", tmp_dir, path_target_on_gcs)
        path_target_on_gcs.fs.put(tmp_dir, path_target_on_gcs, recursive=True)


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
            "--mixture inputbulkrnaseq.txt",
            "--sigmatrix sigmatrix.txt",
            "--sourceGEPs sourcegeps.txt",
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
