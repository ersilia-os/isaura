import csv
from isaura.register.local import IsauraWriter

access_level = "private"

with IsauraWriter(
  input_csv="data.csv",
  model_id="eos3b5e",
  model_version="v1",
  bucket="isaura-public",
  minio_endpoint="http://127.0.0.1:9000",
  minio_access_key="minioadmin",
  minio_secret_key="minioadmin",
  store_directory=".",
  max_rows_per_file=100000,
  checkpoint_every=50000,
  bloom_filename="bloom.pkl",
) as w:
  w.write()
