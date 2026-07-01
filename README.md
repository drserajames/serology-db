# serology-db — local DuckDB / Parquet prototype

A tidy, queryable store of WHO CC influenza serology (HI/FRA/MN/…) built from the
`hidb5` titer databases. Proof-of-concept for ad-hoc analytics that `hidb5`'s
narrow lookup API doesn't expose — and a forward-compatible seed for a future
Postgres-backed upload/query/download service.

> **Code vs data split.** This directory holds **only code** (build scripts,
> `refresh.sh`, README). All output — CSVs, Parquet, and the ~240 MB
> `serology.duckdb` — is written to **`$SERO_OUT`, default
> `acmacs-data/serology-db/`**, which is `.gitignore`d there (it's a large,
> regenerable, WHO-derived artifact). Override the location with the `SERO_OUT`
> env var. Paths written `out/…` below are relative to `$SERO_OUT`.
>
> **Portability:** no paths are hardcoded. Defaults assume `acmacs-data/` and
> `ae/` are siblings of this repo (resolved relative to the scripts); override
> with `ACMACS_DATA`, `AE_BUILD`, and `SERO_OUT` if your layout differs.

## Data lineage

```
acmacs-data/hidb5.{h1,h3,b}.json.xz   (normalised antigen/serum identity + raw titer matrices)
   └─ build_db.py        → out/csv/{antigen,serum,titer_table,titer}.csv   tidy long form
seqdb-{h1,h3,b}.v4.json.xz (HA sequences)
   └─ build_sequences.py → out/csv/match.csv   ag_id→seq_id + aa  (parallel; ~4 min; via ae_backend)
clades.json
   └─ build_clades.py    → out/csv/clade.csv   seq_id→clade      (~4 s; via ae_backend)
locationdb.json.xz
   └─ build_locations.py → out/csv/location.csv  location→country/continent/lat-long (~1 s; pure Python)
        └─ load_duckdb.py → out/serology.duckdb    typed tables + views + indexes
                          → out/parquet/...         portable columnar export
```

Sequences are split into **two stages on different cadences**: the slow *match*
(`ag_id → seq_id`, depends on hidb5+seqdb, ~biweekly) and the near-instant
*clades* (`seq_id → clade`, depends on `clades.json`). `load_duckdb.py` joins
them on `seq_id` into the `sequence` table. See "Updating" for how `refresh.sh`
runs only the stale stage.

We build **from `hidb5`, not the raw xlsx/torg**, because hidb5 already solved the
hard part: antigen/serum identity resolution and passage normalisation. Likewise
sequences **reuse ae's canonical matcher** rather than a hand-rolled name join —
see below.

## Schema (star-ish: fact + dimensions, plus a flat view)

| table | grain | key columns |
|-------|-------|-------------|
| `antigen` | one test antigen | `ag_id`, name, virus_type, lineage, location, year, passage, collection_date |
| `serum` | one antiserum | `sr_id`, serum_id, name, species, lineage, passage |
| `titer_table` | one assay table | `tab_id`, subtype, assay, rbc, lab, virus, table_date |
| `titer` | one antigen×serum reading | `tab_id`, `ag_id`, `sr_id`, titer_raw, titer_kind, titer_value, log_titer |
| `sequence` | one antigen's HA seq | `ag_id`, seq_id, clade, clade_path, aa_length, nuc_length, aa |
| `location` | one resolved place | `location` (raw hidb `O`), canonical, country, continent, division, latitude, longitude |
| `titer_flat` *(view)* | denormalised join | titer + antigen + serum + `antigen_clade` + `antigen_country`/`antigen_continent` |
| `antigen_sequence` *(view)* | antigen ⨝ seq ⨝ geo | antigen columns + seq_id, clade, aa + country/continent/lat/long |

- IDs are `subtype:index` (e.g. `h3:1042`).
- `titer_kind`: `num` (exact), `lt` (`<N`, left-censored), `gt` (`>N`), `other`.
- `log_titer` = log2(value / 10) (acmacs convention). Missing `*` cells are omitted.

Current volume: **268,729 antigens · 5,685 sera · 9,365 tables · 3,581,094 titers
· 102,050 antigen sequences**.

## Sequences & clades — reusing ae's matcher

`build_sequences.py` does **not** re-implement strain-name matching. It drives
`ae_backend`'s seqdb exactly as `ae/py/ae/report/geographic.py` does:

```python
seqdb = ae_backend.seqdb.for_subtype("A(H3N2)")
ref   = seqdb.select_all().filter_name(name=..., reassortant=..., passage=...)[0]
ref.aa; ref.nuc; ref.clades          # after find_clades(clades.json)
```

That C++ matcher applies location-abbreviation expansion (LOCDB_V2), reassortant
and passage handling. Output is keyed by `ag_id`, so it joins to `antigen`/`titer`
exactly (no fuzzy join). Match by name only (passage-agnostic) by default;
`--with-passage` for passage-specific HA.

`filter_name` is CPU-bound, does **not** release the GIL, and names are ~94%
distinct, so `build_sequences.py` fans the work across processes
(`~9.5 min → ~4 min`). It's memory-bandwidth-bound (each match scans the whole
seqdb), so more cores give diminishing returns.

**Incremental cache.** The match is cached by the *stable* natural key
`(subtype, name, reassortant, passage)` in `out/csv/match_cache.csv` (hits **and**
misses), so a refresh only matches what's new or newly matchable:

- cached **hit** → reused (a name hit stays valid; seqdb doesn't drop sequences)
- **new** key → matched
- cached **miss** → re-checked only when seqdb changed, and (by default) only for
  *recent* strains (`year ≥ max_year − 3`) — of 154,725 misses, 131,913 are old
  strains that will never be newly sequenced, so re-trying them is pure waste.

Result: **~4 min first build, then ~5 s** (seqdb unchanged, new keys only) to
**~30 s** (seqdb grew). Flags: `--rematch-all`, `--rematch-misses`,
`--recheck-years N`. Titers and clades stay full-rebuild by design — both are
already ~10 s, so a cache there would add risk for no gain.

**Coverage is ~34–44% of antigens** (h1 44%, h3 37%, b 34%). This is mostly a
genuine fact — many HI antigens were never sequenced — confirmed by the canonical
matcher landing at the same rate a normalised-name join does. Of matched
antigens ~85% get a clade; the rest match a sequence that `clades.json` leaves
unclassified (NULL clade). Most-specific clade in `clade`, full path in `clade_path`.

Needs the built ae extension + reference data (auto-defaulted in the scripts):
`PYTHONPATH`→`ae/build`, `SEQDB_V4`, `LOCDB_V2`, `AC_CLADES_JSON_V2`. The step is
**optional** — if `ae_backend` can't load, the titer DB is unaffected.

> **lab_id fallback — investigated, no gain (2026-07).** ae's chart-level
> `populate_from_seqdb` adds a `select_by_lab_id` fallback when name match fails.
> hidb5 antigens do carry lab_ids (79–87%, `CDC#…` format), **but the local
> seqdb v4 snapshot stores lab *names* with empty id lists** (0 of 214,591 H3
> seqs have a non-empty lab_id). Estimated recovery against this seqdb: **0
> antigens.** It would only help against a seqdb whose lab_ids are populated
> (server-side pipeline). So ~38% is effectively the ceiling here. (Separately,
> `populate_from_seqdb` uses the name-indexed `select_by_name` before
> `filter_name` — a real speed lever, but neither it nor `select_by_lab_id` is
> bound to Python, so both need a `ae_backend` rebuild to use.)

## Locations — resolved via locationdb

`build_locations.py` resolves the raw hidb5 `O` value (our `location` column — a
mix of full names and CDC codes like `WI`, `AG`) to **country, continent,
division, latitude, longitude** using `locationdb.json.xz`. locationdb is
self-contained, so this reads it **directly (no ae dependency)**, following ae's
own resolution order:

```
replacements (spelling) → names (alias) → locations {name: [lat,long,country,division]}
cdc_abbreviations {code: name}   ← fallback for bare CDC 2-letter codes (the ~13% that names miss)
countries {country: idx} → continents[]
```

**~100% coverage:** 99.7% of distinct locations, **99.98% of antigen+serum rows**
(3,838/3,851 locations; the ~13 unresolved are genuine junk like
`cdc-name-without-location`). `location` joins to `antigen`/`serum` on the raw
`location` value; `titer_flat` exposes `antigen_country`/`antigen_continent`, so
you get geographic titer queries (GMT by continent/country, drift by region) and
a bridge to ae's `geo-draw`.

## Build / use

```bash
python3 build_db.py         # hidb5 → out/csv/          (~5 s)
python3 build_locations.py  # locationdb → location.csv (~1 s, pure Python)
python3 build_sequences.py  # seqdb → out/csv/match.csv (~4 min first run, then ~5-30 s; needs ae)
python3 build_clades.py     # clades.json → clade.csv   (~4 s, needs ae)
python3 load_duckdb.py      # CSV → DuckDB + Parquet    (~7 s)
python3 demo_queries.py     # example analytical queries
```

But normally just run `./refresh.sh` — it runs only the stale stages (below).

### Updating (`refresh.sh` — runs only the stale stage)

New tables enter upstream (xlsx → whocc-tables → `.ace` → server-side
`whocc-hidb5-update`) and reach this machine as refreshed
`acmacs-data/hidb5.*.json.xz` + `seqdb-*.v4.json.xz` snapshots. `refresh.sh`
hashes two input groups separately and rebuilds only what each drives:

| Input group | Hash | Drives | Cost | Typical cadence |
|---|---|---|---|---|
| hidb5 + seqdb + locationdb | `DATA` | `build_db` + `build_locations` + `build_sequences` (incremental match) | ~5–30 s (4 min first build) | every 1–2 weeks |
| `clades.json` | `CLADE` | `build_clades` only | **~10 s** | ad-hoc, more often |

```bash
./refresh.sh           # rebuild whatever changed (data and/or clades)
./refresh.sh --force   # rebuild everything
./refresh.sh --no-seq  # titers only; skip the ~4 min sequence match
./refresh.sh --pull    # first pull fresh hidb5 from upstream (needs WHO CC SSH)
```

So a **`clades.json` edit re-derives clades in ~10 s** (build_clades + reload)
without re-running the 4-min match — that was the point of splitting the stages.
Sequence stages are skipped gracefully if `ae_backend` is unavailable.

**Titers: full rebuild (by design).** `build_db`+load is ~12 s and idempotent,
and it's the *safe* choice — hidb5's positional indices (`ag_id`/`tab_id`) are not
stable across regenerations, so an incremental append keyed on them would corrupt.
No cache there buys < 12 s at real risk.

**Match: incremental** — because it's the only expensive stage (~4 min) and it can
key on a *stable natural* key `(subtype, name, reassortant, passage)` rather than
positional indices. See "Sequences & clades" for the cache. State lives in
`out/.manifest` (`DATA`/`CLADE`/`SEQ`) plus `out/csv/match_cache.csv` +
`out/csv/.match_meta`. The same natural-key idea (extended with
lab+assay+rbc+date+content hash for tables) is the prerequisite for a fully
incremental store — needed anyway for the future Postgres upload path.

Query directly with no DB load, language-agnostic, straight off Parquet:

```sql
SELECT subtype, titer_kind, count(*)
FROM read_parquet('../acmacs-data/serology-db/parquet/titer/*/*.parquet', hive_partitioning=true)
GROUP BY ALL;
```

## Known caveats (prototype fidelity)

- **Clipped cells:** ~32.7k titer cells (~0.9%) sit in tables with more titer
  rows/cols than registered hidb indices (antigens/sera absent from the identity
  DB). They are dropped and counted (`clipped_cells` in `build_db.py` output),
  not silently lost. A production build would reconcile these against the `.ace`
  charts.
- **CJK names:** hidb5 uses a non-standard `\U####` escape for Chinese location
  names; `load_json_xz()` normalises it to standard JSON `\u####`.
- **Sequence coverage ~38%:** only antigens with an HA sequence in seqdb get a
  `sequence` row (LEFT-joined, so unsequenced antigens still appear with NULL
  clade). Matching is passage-agnostic by default — egg vs cell HA can differ, so
  use `--with-passage` when that matters.
- **hidb-derived:** titers come from hidb5's merged view, not a fresh parse of
  each `.ace`. Good for analytics; for canonical per-table provenance, go to the
  source charts.

## Path to the eventual service

This schema is deliberately portable: the tidy `titer` fact + dimension tables
lift directly into Postgres as the system-of-record, with this DuckDB/Parquet
layer retained as the fast query + bulk-download tier. The reusable asset across
both is the `ae`/`hidb5` identity normalisation — keep it separable from storage.
