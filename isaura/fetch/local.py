import duckdb, os

MINIO_ENDPOINT = "http://localhost:9000"
BUCKET = "isaura-public"
KEY = "eos3b5e/v1.0.0/tranches_0/eos3b5e_output.csv"
# KEY = "eos5axz/v1/tranches_0/eos5axz_output.csv"
ACCESS_KEY = "minioadmin"
SECRET_KEY = "minioadmin"
print(os.cpu_count())
con = duckdb.connect()
con.execute(f"PRAGMA threads={os.cpu_count() or 8};")
con.execute("INSTALL httpfs; LOAD httpfs;")
con.execute("INSTALL json; LOAD json;")


con.execute(f"""
    SET s3_endpoint='{MINIO_ENDPOINT.replace("http://", "").replace("https://", "")}';
    SET s3_use_ssl={"true" if MINIO_ENDPOINT.startswith("https://") else "false"};
    SET s3_url_style='path';
    SET s3_region='us-east-1';
    SET s3_access_key_id='{ACCESS_KEY}';
    SET s3_secret_access_key='{SECRET_KEY}';
""")


df_exploded = con.sql(f"""
    WITH data AS (
      SELECT * FROM read_csv_auto('s3://{BUCKET}/{KEY}')
    ),
    exploded AS (
      SELECT
        d.*,
      FROM data d
    )
    SELECT *
    FROM exploded
""").fetchdf()
print(df_exploded)


# df = con.sql(f"""
#     SELECT *
#     FROM read_csv_auto(
#         's3://{BUCKET}/{KEY}',
#         HEADER=true,
#         PARALLEL=true
#         -- leave autodetect on since you don't want to list 2000 columns
#     )
#     -- Optional: push down filters to cut work/IO
#     -- WHERE smiles IS NOT NULL
# """).fetchdf()
