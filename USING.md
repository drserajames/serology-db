# Using the serology database

A practical guide to **querying** `serology.duckdb` — the local store of WHO CC
influenza serology (HI / FRA / MN / neutralisation titers, joined to HA sequences,
clades and geography). For how the database is *built* and the design rationale, see
[`README.md`](README.md); this file is just about getting answers out of it.

> **In one sentence:** it's a single [DuckDB](https://duckdb.org) file (plus a Parquet
> mirror) holding ~3.6 M titer measurements across ~9.4 k assay tables, with every
> antigen/serum enriched by strain metadata, sequence/clade and resolved location — so
> most questions are one SQL query away.

---

## 1. Get the file and connect

The database is **not in this repo** (it's WHO-derived data, gitignored). It lives at
`$SERO_OUT/serology.duckdb` — default `~/AC/eu/acmacs-data/serology-db/serology.duckdb`.
Either point at an existing copy, or build one with `./refresh.sh` (see the README; needs
the reference `hidb5`/`seqdb` data). Treat it as a **read-only snapshot** — it's rebuilt
from source each refresh, so don't hand-edit it (changes are overwritten).

Open it read-only (avoids a write lock, lets several readers share the file):

**DuckDB CLI**
```sh
duckdb -readonly ~/AC/eu/acmacs-data/serology-db/serology.duckdb
D SELECT count(*) FROM titer;
```

**Python**
```python
import duckdb
con = duckdb.connect("serology.duckdb", read_only=True)
con.sql("SELECT count(*) FROM titer").show()
df = con.sql("SELECT * FROM titer_flat WHERE subtype = 'h3' LIMIT 1000").df()  # -> pandas
```

**R**
```r
library(duckdb)
con <- dbConnect(duckdb(), "serology.duckdb", read_only = TRUE)
dbGetQuery(con, "SELECT subtype, count(*) FROM titer_table GROUP BY 1")
```

No DuckDB? The same data is exported to Parquet under `$SERO_OUT/parquet/` and readable by
pandas / polars / R-arrow directly (see §6).

A runnable tour of representative queries ships as
[`demo_queries.py`](demo_queries.py) — `python3 demo_queries.py`.

---

## 2. The data model

A star schema: **`titer` is the fact table**, one row per measured antigen×serum cell;
everything else describes its antigen, serum or assay table.

```
             titer_table (assay/lab/date/rbc)         location (country/continent/lat-long)
                    │ tab_id                                 │ location
                    ▼                                        ▼
   antigen ──ag_id──►  titer  ◄──sr_id── serum        antigen.location ─┘
      │                (~3.6M)                         serum.location ───┘
      └─ag_id─► sequence (seq_id, clade, aa)
```

| Object | Grain | Rows | What it holds |
|---|---|---|---|
| `titer` | antigen × serum × table | ~3.6 M | the measurement: `titer_raw`, `titer_kind`, `titer_value`, `log_titer`, `log_titer_thresholded`, `source` |
| `antigen` | antigen | ~270 k | strain name, subtype, lineage, passage, reassortant, location, year, `source` |
| `serum` | serum | ~5.7 k | serum strain + `serum_id`, `species`, subtype, location |
| `titer_table` | assay table | ~9.4 k | `assay`, `rbc`, `lab`, `virus`, `table_date`, subtype |
| `sequence` | antigen | ~102 k | `seq_id`, `clade`, `clade_path`, `aa` (amino-acid string), lengths — only for matched antigens (~38%) |
| `location` | raw location string | ~3.9 k | `country`, `continent`, `division`, `latitude`, `longitude` |

**Two convenience views do the joins for you:**

- **`titer_flat`** — the "one big table": every titer with its assay context, antigen
  (name/lineage/passage/year/clade/country/continent), serum (name/id/species) and both
  log columns. **Start here** — most analysis needs nothing else.
- **`antigen_sequence`** — one row per antigen with its sequence, clade and resolved
  geography attached.

See the exact columns any time with `DESCRIBE titer_flat;` (or `DESCRIBE antigen;`).

---

## 3. Conventions you must know before averaging

These are the traps that turn a valid query into wrong science. Read once.

### Titers are encoded, and some are censored
`titer_raw` is a string. `titer_kind` tells you how to read it:

| `titer_kind` | Example `titer_raw` | Meaning | ~share |
|---|---|---|---|
| `num` | `640` | an exact reading | ~90% |
| `lt` | `<10` | left-censored: below the detection threshold | ~9% |
| `gt` | `>1280` | right-censored: above the top dilution | ~0.5% |

(Truly missing cells are simply **absent** — there is no row, and no `*`. So counting rows
already excludes missing data.)

### Use `log_titer_thresholded` for GMT / SD — not `log_titer`, not `titer_value`
Two log columns exist, and picking the wrong one biases every average:

- `log_titer` = `log2(value / 10)` — the **face value** (a `<10` counts as `10`). This is
  ae's `Titer::logged()`. Fine for display; **biased for averaging**.
- `log_titer_thresholded` = ae's `logged_with_thresholded()` — a `<N` counts as `N/2`
  (`log − 1`), a `>N` as `N×2` (`log + 1`). This is the **canonical convention for
  summary stats**; use it whenever you average.

**Geometric mean titer** is therefore:
```sql
round(pow(2, avg(log_titer_thresholded)) * 10, 1) AS gmt
```
…over `titer_kind IN ('num','lt','gt')` (include the censored readings — that's the point).
Dropping censored titers (`WHERE titer_kind = 'num'`) biases the GMT **upward** because it
discards the low-reactivity tail. Query §5.6 shows the size of that bias.

### Provenance: `source` / `antigen_source`
Rows are tagged `'hidb'` (from the hidb5 merged tables, ~99%) or `'ace_recovered'`
(~32 k cells recovered from source `.ace` charts that hidb5 had dropped). They're
scientifically equivalent; filter `WHERE source = 'hidb'` only if you want a strict
hidb-only view.

### Joins: use the natural keys
`ag_id`, `sr_id`, `tab_id` are **content-based hashes** — stable across rebuilds, so
they're the keys to join on and to reference from other stores. `hidb_id` is the old
positional id, kept only for provenance; don't join on it.

### Missing sequence / location is normal
`sequence` and `location` are **LEFT**-joined in the views: only ~38% of antigens have a
matched sequence (many strains are simply unsequenced), and a few locations don't resolve.
So `antigen_clade`, `antigen_country`, `seq_id` etc. can be `NULL` — filter with
`... IS NOT NULL` when a query needs them.

### Vocabulary
- `subtype`: `'h1'`, `'h3'`, `'b'`.
- `assay`: `'HI'`, `'FRA'`, `'MN'`, `'PRN'`, `'PN'`, `'HINT'`, `'NEUTRALISATION'`, … —
  HI dominates. Passage/RBC context lives in `titer_table`.

---

## 4. Finding things

Strain names aren't shown here (public repo), so these use **bind parameters** — pass the
value from your client instead of pasting a name into the SQL.

```sql
-- all titers for a given antigen strain (bind :name), newest first
SELECT table_date, lab, assay, titer_raw, log_titer_thresholded
FROM titer_flat
WHERE antigen = ?                      -- bind a strain name
ORDER BY table_date DESC;
```
```python
con.execute("SELECT table_date, lab, assay, titer_raw FROM titer_flat "
            "WHERE antigen = ? ORDER BY table_date DESC", ["<strain name>"]).df()
```

Don't know the exact name? Search by fragment (indexes back `antigen.name`):
```sql
SELECT DISTINCT name, subtype, year FROM antigen
WHERE name ILIKE ? ORDER BY year DESC;    -- bind '%fragment%'
```

---

## 5. Query cookbook

All copy-paste-runnable (they pick their own example rows at runtime).

**5.1 Coverage — how much data, by subtype × assay**
```sql
SELECT subtype, assay,
       count(DISTINCT tab_id) AS tables,
       count(*)               AS titers
FROM titer_flat
GROUP BY 1, 2 ORDER BY titers DESC LIMIT 10;
```

**5.2 Reactivity time-series for the most-tested H3 antigen** (runtime-chosen, no name needed)
```sql
SELECT table_date, lab, assay, titer_raw
FROM titer_flat
WHERE subtype = 'h3' AND titer_kind = 'num'
  AND antigen = (SELECT antigen FROM titer_flat
                 WHERE subtype = 'h3' AND titer_kind = 'num'
                 GROUP BY antigen ORDER BY count(*) DESC LIMIT 1)
ORDER BY table_date;
```

**5.3 GMT per serum for the most recent H3 table** (censoring-correct)
```sql
SELECT serum, serum_species, count(*) AS n,
       round(pow(2, avg(log_titer_thresholded)) * 10, 1) AS gmt
FROM titer_flat
WHERE tab_id = (SELECT tab_id FROM titer_table
                WHERE subtype = 'h3' ORDER BY table_date DESC LIMIT 1)
  AND titer_kind IN ('num', 'lt', 'gt')
GROUP BY 1, 2 ORDER BY gmt DESC;
```

**5.4 GMT by antigenic clade** (titers ⨝ seqdb)
```sql
SELECT antigen_clade AS clade,
       count(DISTINCT antigen) AS antigens, count(*) AS titers,
       round(pow(2, avg(log_titer_thresholded)) * 10, 0) AS gmt
FROM titer_flat
WHERE subtype = 'h3' AND antigen_clade IS NOT NULL
  AND titer_kind IN ('num', 'lt', 'gt')
GROUP BY 1 HAVING count(*) > 5000 ORDER BY titers DESC;
```

**5.5 GMT by antigen continent** (titers ⨝ locationdb)
```sql
SELECT antigen_continent AS continent,
       count(DISTINCT antigen) AS antigens, count(*) AS titers,
       round(pow(2, avg(log_titer_thresholded)) * 10, 0) AS gmt
FROM titer_flat
WHERE subtype = 'h3' AND antigen_continent IS NOT NULL
  AND titer_kind IN ('num', 'lt', 'gt')
GROUP BY 1 ORDER BY titers DESC;
```

**5.6 Why the censoring column matters — same H3 HI data, two ways**
```sql
SELECT
  round(pow(2, avg(log_titer)          FILTER (WHERE titer_kind = 'num')) * 10, 1)
    AS gmt_dropping_censored,   -- biased upward
  round(pow(2, avg(log_titer_thresholded) FILTER (WHERE titer_kind IN ('num','lt','gt'))) * 10, 1)
    AS gmt_canonical            -- correct
FROM titer_flat WHERE subtype = 'h3' AND assay = 'HI';
```

**5.7 Provenance — hidb vs cells recovered from source charts**
```sql
SELECT source, count(*) AS titers, count(DISTINCT ag_id) AS antigens
FROM titer GROUP BY 1;
```

**5.8 Censoring / threshold breakdown per subtype** (data-quality view)
```sql
SELECT tb.subtype, t.titer_kind, count(*) AS n,
       round(100.0 * count(*) / sum(count(*)) OVER (PARTITION BY tb.subtype), 1) AS pct
FROM titer t JOIN titer_table tb USING (tab_id)
GROUP BY 1, 2 ORDER BY 1, n DESC;
```

---

## 6. Using the Parquet export (without the DuckDB file)

`$SERO_OUT/parquet/` mirrors the tables; `titer/` is partitioned by subtype.

```python
import polars as pl
ag  = pl.read_parquet("parquet/antigen.parquet")
h3  = pl.read_parquet("parquet/titer/subtype=h3/*.parquet")   # partition pruning
```
```python
# ...or query the Parquet with DuckDB's engine, no .duckdb file needed:
import duckdb
duckdb.sql("SELECT count(*) FROM 'parquet/titer/**/*.parquet'").show()
```

---

## 7. Freshness

The database is a **rebuilt snapshot** of the WHO reference data as of the last
`./refresh.sh`. Row counts drift as sources update. To rebuild or update it (and for the
`--fast-seq` / other options), see the README's *Updating* section. Known fidelity
caveats are listed under the README's *Known caveats*.
```
