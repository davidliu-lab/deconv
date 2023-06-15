import logging

import pandas as pd
from upath.implementations.cloud import GCSPath

import helpers
from helpers.bayesprism import run_and_save_results
from helpers.useful_small_things import make_a_nice_timestamp_of_now

helpers.logging.configure_logging()
logging.getLogger("helpers").setLevel("DEBUG")
logger = logging.getLogger(__name__)
logger.setLevel("DEBUG")

logger.info("hello world")

data_path = GCSPath("gs://liulab/bayesprism")
sc_rnaseq = data_path / "sc_rnaseq.parquet"
sc_rnaseq_annotations = data_path / "annotations.parquet"
bulk_rnaseq = data_path / "bulk_rnaseq.parquet"

if not all([sc_rnaseq.exists(), sc_rnaseq_annotations.exists()]):
    logger.debug("loading jerby-arnon data from GCS")
    df_sc_rnaseq, df_sc_metadata = helpers.datasets.jerby_arnon.load_jerby_arnon_hg19_tpm(
        subset=True
    )
    logger.debug("saving data to GCS")
    df_sc_rnaseq.to_parquet(sc_rnaseq)
    df_sc_metadata.to_parquet(sc_rnaseq_annotations)

if not bulk_rnaseq.exists():
    logger.debug("loading bulk rnaseq data from GCS")
    run_path = GCSPath(
        "gs://liulab/differential_composition_and_expression/20230505_21h41m44s/experiment_id=000/malignant_means=0.55,0.85/log2_fc=-1.50/run_id=00/"
    )
    bulk_rnaseq_list = [
        pd.read_parquet(run_path / group_id / "bulk_rnaseq.parquet") for group_id in ["a", "b"]
    ]
    df_bulk_rnaseq_all = pd.concat(bulk_rnaseq_list, axis=1)
    logger.debug("saving data to GCS")
    df_bulk_rnaseq_all.to_parquet(bulk_rnaseq)

target = GCSPath("gs://liulab/bayesprism") / make_a_nice_timestamp_of_now()

run_and_save_results(sc_rnaseq, sc_rnaseq_annotations, bulk_rnaseq, target)

# execute shell command "ls -la" and print result
# command = f"gsutil ls -R {target}"
# print(subprocess.check_output(command, shell=True).decode("utf-8"))
