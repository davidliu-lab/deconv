# using duckdb

https://duckdb.org/docs/guides/python/sql_on_arrow

https://arrow.apache.org/docs/python/dataset.html

```python
# link to parquet files using an Arrow Dataset
my_arrow_dataset = ds.dataset(str(base_path / 'parquet_folder/'))

# query the Apache Arrow Dataset "my_arrow_dataset" and return as an Arrow Table
results = con.execute("SELECT * FROM my_arrow_dataset WHERE i = 2").arrow()
```
