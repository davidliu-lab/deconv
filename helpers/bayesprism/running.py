import json
import logging
import tempfile

import docker
import pandas as pd
import upath
from upath.implementations.cloud import GCSPath

from .formatting import (
    format_bulk_rnaseq,
    format_sc_cell_states,
    format_sc_cell_types,
    format_sc_rnaseq,
)

logger = logging.getLogger(__name__)


def run(
    df_sc_rnaseq: pd.DataFrame,
    df_sc_rnaseq_annotations: pd.DataFrame,
    df_bulk_rnaseq: pd.DataFrame,
    target_path: GCSPath,
):
    # create tempdir
    # format and save input data to local parquet files in tempdir/input_files
    # run bayesprism in tempdir
    #   - bayesprism saves outputs in tempdir
    # reformat and save all data to target_path
    # this approach...
    #   - avoids auth issues with docker
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = upath.UPath(tmp_dir)
        format_sc_rnaseq(df_sc_rnaseq, df_sc_rnaseq_annotations).write_parquet(
            tmp_path / "sc_rnaseq.parquet"
        )
        format_sc_cell_types(df_sc_rnaseq_annotations).write_parquet(
            tmp_path / "sc_rnaseq_cell_types.parquet"
        )
        format_sc_cell_states(df_sc_rnaseq_annotations).write_parquet(
            tmp_path / "sc_rnaseq_cell_states.parquet"
        )
        format_bulk_rnaseq(df_bulk_rnaseq).write_parquet(tmp_path / "bulk_rnaseq.parquet")
        run_kwargs = dict(
            auto_remove=True,
            detach=True,
            volumes=[f"{tmp_dir}:/running_bayesprism"],
            working_dir="/running_bayesprism",
        )
        command_arguments = " ".join(
            [
                f"--sc_rnaseq_uri sc_rnaseq.parquet ",
                f"--sc_rnaseq_cell_types_uri sc_rnaseq_cell_types.parquet "
                f"--sc_rnaseq_cell_states_uri sc_rnaseq_cell_states.parquet "
                f"--bulk_rnaseq_uri bulk_rnaseq.parquet ",
                # "--window ",  # bayesprism params?
            ]
        )
        logger.debug("running docker run with %s", command_arguments)
        client = docker.from_env()
        container = client.containers.run(
            "grisaitis/bayesprism:latest", command_arguments, **run_kwargs
        )
        for message in container.logs(follow=True, stream=True):
            print(message.decode("utf-8"), end="")
        container.wait()
        target_path.fs.put(tmp_dir, target_path, recursive=True)


def run_and_save_results(
    sc_rnaseq_path: GCSPath,
    sc_rnaseq_cell_types_path: GCSPath,
    sc_rnaseq_cell_states_path: GCSPath,
    bulk_rnaseq_path: GCSPath,
    target_path: GCSPath,
):
    with tempfile.TemporaryDirectory() as tmp_dir:
        logger.debug("running BayesPrism in %s", tmp_dir)
        run_kwargs = dict(
            auto_remove=True,
            detach=True,
            # user=f"{os.getuid()}:{os.getgid()}",
            volumes=[
                f"{tmp_dir}:/bayesprism/data",
                # "/home/jupyter/deconv/helpers/bayesprism/BayesPrism:/bayesprism/BayesPrism:ro",
            ],
            # working_dir="/bayesprism",  # should be set in Dockerfile
        )
        logger.debug("run_kwargs:\n%s", json.dumps(run_kwargs, indent=2, sort_keys=True))
        command_arguments = " ".join(
            [
                f"--sc_rnaseq_uri {sc_rnaseq_path} ",
                f"--sc_rnaseq_cell_types_uri {sc_rnaseq_cell_types_path} "
                f"--sc_rnaseq_cell_states_uri {sc_rnaseq_cell_states_path} "
                f"--bulk_rnaseq_uri {bulk_rnaseq_path} ",
                # "--window ",  # bayesprism params?
            ]
        )
        logger.debug("running docker run with %s", command_arguments)
        client = docker.from_env()
        container = client.containers.run(
            "grisaitis/bayesprism:latest", command_arguments, **run_kwargs
        )
        for message in container.logs(follow=True, stream=True):
            print(message.decode("utf-8"), end="")
        container.wait()
        target_path.fs.put(tmp_dir, target_path, recursive=True)
