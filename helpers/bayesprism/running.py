import json
import logging
import tempfile

import docker
import gcsfs
import pandas as pd
import upath
from upath.implementations.cloud import GCSPath

logger = logging.getLogger(__name__)


def run_and_save_results(
    reference_sc_rnaseq: GCSPath,
    reference_sc_rnaseq_annotations: GCSPath,
    bulk_rnaseq: GCSPath,
    target_on_gcs: GCSPath,
):
    with tempfile.TemporaryDirectory() as tmp_dir:
        logger.debug("running BayesPrism in %s", tmp_dir)
        run_kwargs = dict(
            auto_remove=True,
            detach=True,
            # user=f"{os.getuid()}:{os.getgid()}",
            volumes=[
                f"{tmp_dir}:/bayesprism/data",
                "/home/jupyter/deconv/helpers/bayesprism/BayesPrism:/bayesprism/BayesPrism:ro",
            ],
            working_dir="/bayesprism",
        )
        logger.debug("run_kwargs:\n%s", json.dumps(run_kwargs, indent=2, sort_keys=True))
        command_arguments = " ".join(
            [
                f"--reference_sc_rnaseq_uri {reference_sc_rnaseq} ",
                f"--reference_sc_rnaseq_annotations_uri {reference_sc_rnaseq_annotations} "
                f"--bulk_rnaseq_uri {bulk_rnaseq} ",
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
        target_on_gcs.fs.put(tmp_dir, target_on_gcs, recursive=True)
