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

# gcs, path = fs.FileSystem.from_uri("gs://liulab")
# print(gcs)
# print(path)

gene_stats = ds.dataset(
    "gs://liulab/differential_composition_and_expression/copied/20230505_21h41m44s/deg_analysis/",
    format="parquet",
    partitioning="hive",
)

relation_gene_stats = duckdb.from_arrow(gene_stats)

# parquet files of gene_stats
# example_path = GCSPath(
#     "gs://liulab/differential_composition_and_expression/20230505_21h41m44s/experiment_id=032/malignant_means=0.65,0.75/log2_fc=0.50/run_id=00/deg_analysis/bulk/gene_stats.parquet"
# )

# base_path = GCSPath("gs://liulab/differential_composition_and_expression/20230505_21h41m44s")
# files = list(base_path.glob("**/deg_analysis/malignant_cibersortx/gene_stats.parquet"))
# file_keys = list(map(lambda p: p.path[1:], files))

# import dask.dataframe as dd
# ddf_malignant_cibersortx_gene_stats = dd.read_parquet(files, engine="pyarrow")
# ddf_malignant_cibersortx_gene_stats
