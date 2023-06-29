import concurrent.futures
import logging
import re

from google.cloud import storage

logger = logging.getLogger(__name__)


def rename_files(bucket_name, prefix):
    storage_client = storage.Client()
    bucket = storage_client.get_bucket(bucket_name)

    # List all objects with the given prefix
    blobs = bucket.list_blobs(prefix=prefix)

    for blob in blobs:
        # Extract the original name
        original_name = blob.name
        # remove the 'experiment_id=.*/' part from the name, using re.sub, escaping the / with \
        new_name = re.sub(r"experiment_id=.*/malignant_means", "malignant_means", original_name)

        if new_name != original_name:
            # print(original_name)
            # print(new_name)
            # Rename the blob
            new_blob = bucket.rename_blob(blob, new_name)
            print(f"Renamed {original_name} to {new_blob.name}")


def rename_file(blob):
    # logger.info("blob name: %s", blob.name)
    original_name = blob.name
    bucket = blob.bucket
    new_name = re.sub(r"experiment_id=.*/malignant_means", "malignant_means", original_name)
    if new_name != original_name:
        # return f"Renamed {original_name} to {new_name}"
        new_blob = bucket.rename_blob(blob, new_name)
        return f"Renamed {original_name} to {new_blob.name}"
    else:
        return f"Skipped {original_name}"


def rename_files_2(bucket_name, prefix):
    storage_client = storage.Client()
    bucket = storage_client.get_bucket(bucket_name)

    blobs = bucket.list_blobs(prefix=prefix)

    with concurrent.futures.ThreadPoolExecutor() as executor:
        # Submit rename_file function for each blob
        rename_futures = [executor.submit(rename_file, blob) for blob in blobs]
        logging.info("number of rename_futures: %s", len(rename_futures))

        # Collect results as they become available
        for future in concurrent.futures.as_completed(rename_futures):
            try:
                result = future.result()
                print(result)
            except Exception as e:
                print(f"Error: {e}")


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(threadName)s %(levelname)s %(name)s %(message)s"
    )
    # rename_files("liulab", "differential_composition_and_expression/20230615_01h18m52s/deg_analysis/experiment_id=")
    rename_files_2(
        "liulab",
        "differential_composition_and_expression/20230615_01h18m52s/",
    )


# # Usage example
# bucket_name = "liulab"
# prefix = "differential_composition_and_expression/20230615_01h18m52s/deg_analysis/"
# rename_files(bucket_name, prefix)
