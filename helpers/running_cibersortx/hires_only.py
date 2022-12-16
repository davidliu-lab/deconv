import json
import logging
import pathlib
import tempfile
from typing import Union

import docker
import pandas as pd
import upath

from helpers.running_cibersortx.creating_input_files import (
    create_csx_fractions_tsv,
    create_csx_mixtures_tsv,
)

logger = logging.getLogger(__name__)


def run_and_upload_from_dataframes(
    df_bulk_rnaseq: pd.DataFrame,
    df_cibersort_results: pd.DataFrame,
    path_target_on_gcs: Union[upath.UPath, pathlib.Path],
):
    with tempfile.TemporaryDirectory() as tmp_dir:
        logger.debug("tmp_dir: %s", tmp_dir)
        logger.debug("watch -n 0.1 tree -ghpu %s", tmp_dir)
        csx_path = pathlib.Path(tmp_dir)
        (csx_path / "data").mkdir()
        (csx_path / "outdir").mkdir()
        create_csx_mixtures_tsv(df_bulk_rnaseq, csx_path / "data" / "bulkrnaseq.txt")
        create_csx_fractions_tsv(df_cibersort_results, csx_path / "outdir" / "fractions.txt")
        run(tmp_dir)
        logger.debug("copying tmp_dir %s to %s", tmp_dir, path_target_on_gcs)
        path_target_on_gcs.fs.put(tmp_dir, path_target_on_gcs, recursive=True)


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
    logger.debug("run_kwargs:\n%s", json.dumps(run_kwargs, indent=2, sort_keys=True))
    command_arguments = " ".join(
        [
            "--username lyronctk@stanford.edu",
            "--token dfeba2c8b9d61daebee5fa87026b8e56",
            # "--username william_grisaitis@dfci.harvard.edu",
            # "--token fd214a5e8fae301ca42c195f21ef2602",
            "--mixture bulkrnaseq.txt",  # <file_name>  Mixture matrix [required]
            # "--verbose TRUE",  # not sure this is even an option?
            # "--groundtruth ",   # <file_name>  Ground truth GEPs [same labels as classes] [optional; default: none]            "--mixture bulkrnaseq.tsv",
            "--cibresults fractions.txt",
            # "--degclasses ",  # <file_name>  Run on two classes, specified by 1, 2, 0=skip [default: none]
            # "--filtered ", # <file_name>  Filtered GEPs from CIBERSORTxGEP [default: none]
            # "--window ",  # <int>   Window size for deconvolution [default: No. of cell types x 4]
            # "--useadjustedmixtures ",  # <bool>  If doing B-mode batch correction, use adjusted mixtures for GEP imputation [default: FALSE]
        ]
    )
    logger.debug("docker run --rm cibersortx/hires:latest %s", command_arguments)
    client = docker.from_env()
    container = client.containers.run("cibersortx/hires:latest", command_arguments, **run_kwargs)
    for message in container.logs(follow=True, stream=True):
        print(message.decode("utf-8"), end="")
    container.wait()
