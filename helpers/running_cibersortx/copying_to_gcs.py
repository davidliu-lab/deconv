import logging
import pathlib
import shutil

from upath import UPath

logger = logging.getLogger(__name__)


def copy_file_maybe_in_the_cloud_to_local_path(source_uri: str, target_path: pathlib.Path):
    path = UPath(source_uri)
    try:
        path.fs.get(path, target_path)
    except AttributeError:
        shutil.copy(path.resolve(), target_path.resolve())
