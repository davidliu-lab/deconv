import dask.dataframe as dd
import pandas as pd
import numpy as np
import dask.distributed


if __name__ == "__main__":
    rng = np.random.default_rng(42)
    client = dask.distributed.Client()
    print(client)
    df = pd.DataFrame(
        {
            "a": rng.integers(0, 100, 1000),
            "b": rng.integers(0, 100, 1000),
        }
    )
    ddf = dd.from_pandas(df, npartitions=4)
    print(ddf.head())

    ddf = ddf.assign(c=lambda x: x.a + x.b)
    print(ddf.head())

    ddf = ddf.assign(d=lambda x: x.c * x.b)
    print(ddf.head())

    print(ddf.compute())

    client.close()
