import logging
from concurrent.futures import ThreadPoolExecutor

from upath.implementations.cloud import GCSPath

logger = logging.getLogger(__name__)

base_path = GCSPath("gs://liulab/differential_composition_and_expression/20230505_21h41m44s")
# files = list(base_path.glob("**/deg_analysis/malignant_cibersortx/gene_stats.parquet"))
files = list(base_path.glob("**/gene_stats.parquet"))

new_base = GCSPath("gs://liulab/differential_composition_and_expression/copied")

thread_pool = ThreadPoolExecutor(max_workers=10)

from google.cloud import storage

bucket_liulab = storage.Client().bucket("liulab")


def copy_file(origin: GCSPath, destination: GCSPath):
    assert type(origin) == GCSPath
    assert type(destination) == GCSPath
    source_blob = bucket_liulab.blob(origin.path[1:])
    assert source_blob.exists(), f"\n gsutil ls -l gs://liulab{origin.path}"
    logger.debug("copying %s to %s", origin.path[1:], destination.path[1:])
    blob_copy = bucket_liulab.copy_blob(source_blob, bucket_liulab, destination.path[1:])
    logger.debug("done copying")
    assert blob_copy.exists()


if __name__ == "__main__":
    logging.basicConfig(
        level="DEBUG",
        format="%(asctime)s - %(name)s - %(threadName)s - %(levelname)s - %(message)s",
    )

    thread_futures = []

    for path in files:
        # print(path)
        # print(path.path)
        new_path = new_base / path.parts[2]  # timestamp
        new_path /= path.parts[-3]  # parquet dataset name
        new_path /= f"origin={path.parts[-2]}"  # origin (bulk or malignant_cibersortx)
        for part in path.parts[-7:-3]:
            new_path /= part
        new_path /= path.parts[-1]  # filename
        # print(new_path)
        copy_file(path, new_path)
    #     thread_future = thread_pool.submit(copy_file, path, new_path)
    #     thread_futures.append(thread_future)
    #     # break

    # for thread_future in thread_futures:
    #     print(thread_future.result(timeout=5))
