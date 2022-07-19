import glob
import os
import pathlib

from google.cloud import storage


def copy_local_directory_to_gcs(local_path: str, bucket: storage.Bucket, gcs_path: str):
    """Recursively copy a directory of files to GCS.

    local_path should be a directory and not have a trailing slash.
    """
    assert os.path.isdir(local_path)
    for local_file in glob.glob(local_path + "/**"):
        if not os.path.isfile(local_file):
            continue
        remote_path = os.path.join(gcs_path, local_file[1 + len(local_path) :])
        blob = bucket.blob(remote_path)
        blob.upload_from_filename(local_file)


def copy_file_maybe_in_the_cloud_to_local_path(
    source_uri: str, target_path: pathlib.Path
):
    path = AnyPath(source_uri)
    try:
        path.copy(target_path)
    except AttributeError:
        shutil.copy(path.resolve(), target_path.resolve())
