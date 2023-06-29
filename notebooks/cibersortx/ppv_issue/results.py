import logging

import duckdb
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq
import upath
from pyarrow import fs
from upath.implementations.cloud import GCSPath

logger = logging.getLogger(__name__)

logger.debug("loading module %s", __name__)


gene_stats = ds.dataset(
    "gs://liulab/differential_composition_and_expression/copied/20230505_21h41m44s/deg_analysis/",
    format="parquet",
    partitioning="hive",
)
