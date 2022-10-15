import dask.dataframe as dd
import pandas as pd
import numpy as np
import dask.distributed


if __name__ == "__main__":
    rng = np.random.default_rng(42)
    client = dask.distributed.Client()
    print(client)
