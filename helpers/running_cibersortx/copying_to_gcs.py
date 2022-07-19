import glob
import logging
import os
import pathlib
import shutil
import time

from cloudpathlib import AnyPath, GSPath
from google.cloud import storage

import helpers

logger = logging.getLogger(__name__)


def copy_local_directory_to_gcs(local_path: str, bucket: storage.Bucket, gcs_path: str):
    path = pathlib.Path(local_path)
    for i in path.glob("**/*"):
        if i.is_dir():
            continue
        t = time.time()
        target_gcs_path = GSPath(gcs_path) / i.relative_to(local_path)
        blob = bucket.blob(target_gcs_path.blob)
        blob.upload_from_filename(i.resolve())
        logger.debug(
            f"Took {time.time() - t:4.3f} seconds to upload {i} to {target_gcs_path}"
        )


def copy_file_maybe_in_the_cloud_to_local_path(
    source_uri: str, target_path: pathlib.Path
):
    path = AnyPath(source_uri)
    try:
        path.copy(target_path)
    except AttributeError:
        shutil.copy(path.resolve(), target_path.resolve())


if __name__ == "__main__":
    logger.setLevel("DEBUG")
    helpers.logging.configure_logging()
    bucket = storage.Client().bucket("liulab")
    local_dir = "/home/jupyter/deconv/tmp"
    copy_local_directory_to_gcs(local_dir, bucket, gcs_path="gs://liulab/tmp/stuff2/")
