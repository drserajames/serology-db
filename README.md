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

## Guardrail — pre-commit hook

This is a **public** repo, so a tracked hook (`hooks/pre-commit`) blocks any
commit whose staged changes look like real assay data — strain names, serum/lab
IDs, GISAID accessions, or amino-acid/nucleotide sequence runs. Activate it once
per clone:

```bash
git config core.hooksPath hooks
```

Bypass a genuine false positive with `git commit --no-verify`. The hook is a
backstop; the primary protection is that all data output lives in the gitignored
`$SERO_OUT`, never in this tree.

## Data lineage

```
acmacs-data/hidb5.{h1,h3,b}.json.xz   (normalised antigen/serum identity + raw titer matrices)
   └─ build_db.py        → out/csv/{antigen,serum,titer_table,titer}.csv   tidy long form
        (ids = content-based natural keys via natural_keys.py; see "Natural keys")
whocc-tables/*.ace  (canonical source charts)
   └─ build_clipped.py   → out/csv/recovered_{antigen,serum,titer}.csv  cells hidb dropped (opt; via ae_backend)
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
| `antigen` | one test antigen | `ag_id`, hidb_id, name, virus_type, lineage, location, year, passage, collection_date, source |
| `serum` | one antiserum | `sr_id`, hidb_id, serum_id, name, species, lineage, passage, source |
| `titer_table` | one assay table | `tab_id`, hidb_id, subtype, assay, rbc, lab, virus, table_date |
| `titer` | one antigen×serum reading | `tab_id`, `ag_id`, `sr_id`, titer_raw, titer_kind, titer_value, log_titer, log_titer_thresholded, source |
| `sequence` | one antigen's HA seq | `ag_id`, seq_id, clade, clade_path, aa_length, nuc_length, aa |
| `location` | one resolved place | `location` (raw hidb `O`), canonical, country, continent, division, latitude, longitude |
| `titer_flat` *(view)* | denormalised join | titer + antigen + serum + `antigen_clade` + `antigen_country`/`antigen_continent` |
| `antigen_sequence` *(view)* | antigen ⨝ seq ⨝ geo | antigen columns + seq_id, clade, aa + country/continent/lat/long |

- IDs are **content-based natural keys** (`natural_keys.py`): `{sub}:a:{hash}`
  (antigen), `{sub}:s:{hash}` (serum), `{sub}:t:{hash}` (table). The antigen/serum
  hash is ae's canonical *designation* — name + reassortant + annotations + passage;
  the table hash is lab + assay + rbc + date + virus + a **reorder-invariant** hash
  of its (antigen, serum, titer) content. Unlike hidb's positional indices these are
  **stable across regenerations** (same identity → same id regardless of array
  order), which is the prerequisite for incremental updates and the Postgres store.
  The old positional id (`{sub}:{i}`) is retained as `hidb_id` for provenance.
  Recovered antigens/sera (see "Clipped cells") use `{sub}:ra:{hash}` / `{sub}:rs:{hash}`
  in the same scheme. See "Natural keys" below.
- `source`: `hidb` (the bulk) or `ace_recovered` (clipped cells recovered from the
  source charts). Filter on it to include/exclude recovered data; `titer_flat`
  exposes both `source` (of the titer) and `antigen_source`.
- `titer_kind`: `num` (exact), `lt` (`<N`, left-censored), `gt` (`>N`), `other`.
- `log_titer` = log2(value / 10) (acmacs convention). Missing `*` cells are omitted.
  This is exactly ae's `Titer::logged()` (verified byte-for-byte against the C++
  toolkit on real charts — see `tests/test_ae_fidelity.py`).
- `log_titer_thresholded` = the **aggregation-correct** log, matching ae's
  `Titer::logged_with_thresholded()`: identical to `log_titer` for exact readings,
  but a left-censored `<N` counts as `N/2` (`log−1`) and a right-censored `>N` as
  `N×2` (`log+1`). **Use this column for GMT/SD/averaging.** See "Censoring" below.

Current volume (hidb): **268,728 antigens · 5,684 sera · 9,365 tables · 3,581,094
titers · 102,050 antigen sequences** (antigens/sera are one lower than the raw hidb
record count because natural-key dedup collapses one byte-identical duplicate each),
plus **+1,021 antigens · +36 sera · +32,043 titers** recovered from source charts
(`source='ace_recovered'`) when `build_clipped` runs → 269,749 antigens · 5,720 sera
· 3,613,137 titers total.

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
> (server-side pipeline). So ~38% is effectively the ceiling here.
>
> **pybind speed lever — evaluated, not worth it (2026-07).** The match's per-key
> cost (~2.4 ms) is dominated by rebuilding the full `select_all()` selection each
> call (~2.0 ms; memory-bandwidth-bound), not the name scan. Exposing ae's
> `select_by_name` to Python would *not* help: it is an **O(n) linear scan**, not
> index-backed (`cc/sequences/seqdb.cc:168` — it ignores the class's own O(log n)
> `find_by_name` binary search). A real speedup would need a C++ *logic* change
> (reimplement `select_by_name` on `find_by_name`) plus a rebuild of the fragile ae
> tree and a fork coupling — for a one-time cold-build saving, since the match is
> already incrementally cached (~5–30 s warm). Cheaper alternative if cold-build
> speed ever matters: a **pure-Python name index** — iterate the seqdb once (~0.15 s
> for 241k H3 seqs via the bound `.name()`), then dict-lookup each antigen. Measured
> **100 % hit/miss agreement** with `filter_name` and 136/137 exact seq_id (the one
> diff is a same-strain ranking tie-break), turning the ~4 min cold match into ~1 s
> with **no ae rebuild**. Left unimplemented (caching already covers the common case).

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
python3 build_clipped.py    # .ace → recovered_*.csv    (~3 s, optional; needs ae + whocc-tables)
python3 build_locations.py  # locationdb → location.csv (~1 s, pure Python)
python3 build_sequences.py  # seqdb → out/csv/match.csv (~4 min first run, then ~5-30 s; needs ae)
python3 build_clades.py     # clades.json → clade.csv   (~4 s, needs ae)
python3 load_duckdb.py      # CSV → DuckDB + Parquet    (~7 s)
python3 demo_queries.py     # example analytical queries
```

But normally just run `./refresh.sh` — it runs only the stale stages (below).

### Testing

```bash
python3 -m pip install -r requirements-dev.txt   # duckdb + pytest
python3 -m pytest                                # from the repo root
```

Three layers, all hermetic (no WHO data, no network):

- **Unit** (`test_parsing.py`, `test_natural_keys.py`) — pins `build_db`'s pure
  reshape functions (`parse_titer` kind/value/log, strain-name reconstruction) and
  the natural-key scheme (determinism, field discrimination, reorder-invariance).
- **Fixture regression** (`test_build_db_fixture.py`) — runs the real
  `build_db.main()` over a tiny synthetic `hidb5.h1.json.xz` in a tmp dir and
  asserts the emitted CSVs exactly, including the clipped-cell tally. This is the
  guard for changes to the reshape/clipping logic.
- **DB invariants** (`test_db_invariants.py`) — read-only checks against the
  built `serology.duckdb`: FK integrity (no orphan titers/sequences), unique PKs,
  titer-kind partition, `log_titer == log2(value/10)`, location resolution ≥ 98%,
  sequence-coverage band, and that the `titer_flat` / `antigen_sequence` views
  preserve row counts. Skips cleanly (not fails) if the DB hasn't been built.

The synthetic fixtures contain only obviously-fake values and live in a tmp dir —
never committed. Run the invariant layer after any `./refresh.sh` to confirm the
refresh didn't silently corrupt the store.

**CI.** `.github/workflows/ci.yml` runs the *hermetic* layers (unit + fixture) on
Python 3.12–3.14 for every push/PR — CI has no WHO data, no built `ae_backend`,
and no charts, so the invariant/fidelity layers skip there (27 run, 18 skip). A
green run means the reshape/parse logic is sound; run the full suite locally.

### Environment & reproducibility

| Piece | Value / source |
|---|---|
| Python | 3.14 (`/opt/homebrew/bin/python3`); code is stdlib + `duckdb` only |
| Runtime deps | `requirements.txt` (`duckdb~=1.5`) — `pip install -r requirements.txt` |
| Dev/test deps | `requirements-dev.txt` (adds `pytest~=9.1`) |
| Sequence/clade deps | `ae_backend` (built from `~/AC/eu/ae`, reached via `PYTHONPATH=<ae>/build`) — **not a pip package**; optional (titer DB builds without it) |
| Env overrides | `ACMACS_DATA`, `AE_BUILD`, `SERO_OUT`, `SERO_CSV_DIR`; sequence stages also read `SEQDB_V4`, `LOCDB_V2`, `AC_CLADES_JSON_V2` (all auto-defaulted relative to the repo) |

Quick preflight before a build:

```bash
python3 -c "import duckdb; print('duckdb', duckdb.__version__)"
PYTHONPATH=<ae>/build python3 -c "import ae_backend; print('ae_backend ok')"  # optional
```

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

## Censoring — how left/right-censored titers are logged

Assay titers are frequently censored: `<10` (below the detection limit) or
`>2560` (off the top of the dilution series). In this data **9.3% of readings are
left-censored** (`lt`) and 0.6% right-censored (`gt`) — not negligible for summary
statistics. There are two logged columns, mirroring ae's two purpose-specific
conventions (validated against the C++ `Titer` class in `tests/test_ae_fidelity.py`):

| column | ae function | `<10` → | `>1280` → | use for |
|---|---|---|---|---|
| `log_titer` | `logged()` | 0 (face value) | 7 | display; regular-titer map distances |
| `log_titer_thresholded` | `logged_with_thresholded()` | −1 (= 5) | 8 (= 2560) | **GMT / SD / averaging** |

The choice is not cosmetic. H3 HI GMT computed three ways over the same readings:

| treatment | GMT |
|---|---|
| drop censored (`titer_kind='num'` only) | **205** |
| face value (`avg(log_titer)`, all kinds) | 182 |
| canonical (`avg(log_titer_thresholded)`) | **174** |

Dropping censored titers biases GMT *upward* by ~18% here, because it discards the
low `<N` readings. The demo queries (and any GMT you write) should use
`log_titer_thresholded` and include `lt`/`gt` — see `demo_queries.py` query 8.
(For *map optimization* ae treats `<N` as a one-sided inequality constraint and
drops `>N`; that lives in the optimizer, not in this analytics store.)

## Known caveats (prototype fidelity)

- **Clipped cells — recovered (`build_clipped.py`).** 32,718 titer cells (~0.9%)
  across 1,388 tables sit in matrix rows/cols beyond the registered hidb indices —
  antigens/sera hidb5 itself could not resolve to its identity DB (ae's own hidb
  reader drops them identically). `build_db` drops and counts them; the optional
  `build_clipped` stage then **recovers them from the canonical source charts** in
  `whocc-tables/`. It maps each clipped table to its `.ace` and — only when the
  registered cells are *byte-identical* to the source (strict alignment, so a
  wrong strain identity is never attached) — re-adds the dropped rows/cols as
  `antigen`/`serum`/`titer` rows tagged `source='ace_recovered'`. Result:
  **32,043 / 32,718 cells (97.9%) recovered** across 1,352 tables (1,021 distinct
  recovered antigens, 36 sera). The residual is 36 tables with no aligned source
  chart on disk (~1.6%) plus 164 true replicate-grain cells dropped to keep the
  `(table, antigen, serum)` grain unique. Recovery is optional (needs `ae_backend`
  + `whocc-tables`); a titers-only build simply omits it.
- **CJK names:** hidb5 uses a non-standard `\U####` escape for Chinese location
  names; `load_json_xz()` normalises it to standard JSON `\u####`.
- **Sequence coverage ~38%:** only antigens with an HA sequence in seqdb get a
  `sequence` row (LEFT-joined, so unsequenced antigens still appear with NULL
  clade). Matching is passage-agnostic by default — egg vs cell HA can differ, so
  use `--with-passage` when that matters.
- **hidb-derived:** titers come from hidb5's merged view, not a fresh parse of
  each `.ace`. Good for analytics; for canonical per-table provenance, go to the
  source charts.

## Natural keys (`natural_keys.py`)

hidb5's antigen/serum/table ids are **positional array indices** — `h3:1042` is
"the 1042nd h3 antigen in *this* build". That index is **not stable across hidb
regenerations**: the array reorders, so the same id can mean a different strain
next fortnight. That's why titers are a full rebuild (an incremental append keyed
on positional ids would silently corrupt) and it blocks any incremental/Postgres
store.

So ids are **derived from identity content** instead (shared by `build_db` and
`build_clipped` via `natural_keys.py`, so the two can't drift):

| entity | key | hash input |
|---|---|---|
| antigen | `{sub}:a:{h}` | name, reassortant, annotations, passage (ae "designation") |
| serum | `{sub}:s:{h}` | serum_id, name, reassortant, annotations, passage, species |
| table | `{sub}:t:{h}` | lab, assay, rbc, date, virus, **reorder-invariant** content hash |

The table content hash is the sorted set of `(ag_key, sr_key, titer)` triples, so a
row/column reshuffle across a regeneration leaves the id unchanged; the `{sub}:`
prefix keeps subtype partitioning working. The old positional id is preserved in
`hidb_id` for provenance/debugging.

**Validated** against the current hidb5 (see `tests/test_natural_keys.py` +
`test_db_invariants.test_titer_grain_unique`): over 3.58M titers the
`(table, antigen, serum)` grain has **0 collisions and 0 value conflicts**; the
keys collapse exactly 1 byte-identical duplicate antigen and 1 serum (correct
dedup). Annotations are part of the antigen/serum key precisely because omitting
them merged distinct egg/cell variants and produced grain conflicts.

## Path to the eventual service

This schema is deliberately portable: the tidy `titer` fact + dimension tables
lift directly into Postgres as the system-of-record, with this DuckDB/Parquet
layer retained as the fast query + bulk-download tier. The reusable asset across
both is the `ae`/`hidb5` identity normalisation, now crystallised as the
**content-based natural keys above** — a stable `(table, antigen, serum)` grain is
exactly what an incremental UPSERT needs. Keep `natural_keys.py` separable from
storage.
