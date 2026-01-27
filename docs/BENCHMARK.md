## Isaura Querying Performance Benchmark

This report summarizes the performance improvements achieved in the Isaura querying engine, comparing the original tranche-based DuckDB implementation with optimized versions leveraging Hive-style data partitioning and join-less querying.

---

### Overview

The original DuckDB query engine used a tranche-based structure (e.g., `tranch_0_0`) and performed queries via view registration and inner joins between in-memory DataFrames and Parquet data. While functional, this approach resulted in heavy I/O operations and slow query performance when handling high-dimensional datasets.

#### Original Structure

```
tranches/
  tranch_0_0/
  tranch_0_1/
  tranch_1_0/
  ...
```

#### Improved Hive-Based Structure

```
eos5axz/v1/tranches/
  row=1/
      col=1/
      col=2/
      ...
  row=2/
      col=1/
      ...
```

The Hive-style directory layout improved data partition pruning and I/O locality, allowing DuckDB to leverage predicate pushdown for faster queries.

---

### Benchmark Setup

| Parameter | Value    |
| --------- | -------- |
| Model     | eos5axz  |
| Data type | int      |
| Dimension | 2048     |
| Data size | 10k      |
| Caching   | Disabled |

#### Baseline Performance (Tranche-Based)

| Configuration | Elapsed (s) | Rate (rows/s) |
| ------------- | ----------- | ------------- |
| Without Hive  | 341.21      | 29.1          |
| With Hive     | 276.08      | 36.2          |

---

### Query Optimization Steps

Further performance gains were achieved by:

* Disabling view registration
* Removing INNER JOIN operations
* Allowing direct `WHERE IN (SELECT * FROM UNNEST(?))` filtering

These changes allowed DuckDB to perform predicate pushdown, reducing Parquet scan size and improving CPU efficiency.

#### Results for eos5axz

| Data Size | Elapsed (s) | Rate (rows/s) |
| --------- | ----------- | ------------- |
| 10k       | 15.97       | 626.1         |
| 30k       | 17.06       | 882.0         |

#### Results for eos8a4x

| Model   | dtype | Dim | Data Size       | Elapsed (s) | Rate (rows/s) |
| ------- | ----- | --- | --------------- | ----------- | ------------- |
| eos8a4x | float | 200 | 10k             | 4.11        | 9722.9        |
| eos8a4x | float | 200 | 10k / 100k scan | 6.29        | 17496.2       |
| eos8a4x | float | 200 | 50k / 100k scan | 7.48        | 14709.3       |

---

### Memory Efficiency

Isaura demonstrated excellent memory efficiency during both read and write operations, largely due to:

* Reduced in-memory intermediate views
* Streaming parquet scans
* Optimized batch sizes

This improvement not only accelerates query time but also enables scalable inference across large molecular datasets.
