"""A package for using deconvolution to study metastatic melanoma"""

__version__ = "0.1"
__author__ = "William Grisaitis"

from . import (
    cell_type_naming,
    columns,
    creating_mixtures,
    datasets,
    logging,
    pseudobulk_evaluation,
)
from .download_ftp_file import download_gz_from_ftp
