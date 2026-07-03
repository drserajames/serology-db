# serology-db — Graduation Roadmap

This repo is a working **prototype** database of WHO CC influenza serology (HI/FRA/MN
titers + HA sequences/clades + geography). The goal captured here is to **graduate it
toward production-grade**: improve data fidelity, add validation, and harden the
engineering — while keeping it working at every step. This is an incremental campaign,
not a rewrite. It doubles as a handoff for anyone (human or agent) picking the work up
cold; start by reading `README.md` end to end (it's the authoritative design doc).

## Where things live

- **Code repo:** `~/AC/eu/serology-db/` — this repo, branch `master`, remote `origin` =
  `github.com/drserajames/serology-db` (**PUBLIC**).
- **Data output:** `$SERO_OUT` (default `~/AC/eu/acmacs-data/serology-db/`) — the
  ~240 MB `serology.duckdb`, CSVs, Parquet. **Gitignored; never committed.**
- **Source reference data:** `~/AC/eu/acmacs-data/` — `hidb5.{h1,h3,b}.json.xz` (merged
  WHO titer tables), `seqdb-{h1,h3,b}.v4.json.xz` (HA sequences), `locationdb.json.xz`,
  `clades.json`.
- **Canonical source charts:** `~/AC/eu/whocc-tables/` (~13.7k `.ace` files) — the ground
  truth for the fidelity work. Read-only.
- **ae toolkit:** `~/AC/eu/ae/` — provides `ae_backend` (pybind11 C++ extension) used for
  sequence matching + clades.

## How to run / verify

- Python is `/opt/homebrew/bin/python3` (3.14); `duckdb` installed via pip
  (`--break-system-packages`). Shell is zsh (no `shopt`).
- **Full refresh:** `cd ~/AC/eu/serology-db && ./refresh.sh` — hash-gated, runs only
  stale stages. First sequence build ~4 min; warm ~seconds.
- **Confirm ae works:**
  `PYTHONPATH=~/AC/eu/ae/build python3 -c "import ae_backend; ae_backend.seqdb.for_subtype('A(H3N2)')"`
  If it fails, the ae build may be stale (it has broken before — a worktree with compile
  errors). Check `git -C ~/AC/eu/ae worktree list` and that
  `ae/build/ae_backend.cpython-314-darwin.so` exists; rebuild with ninja if needed. The
  sequence/clade stages skip gracefully if ae is absent (the titer DB still builds).
- **Sanity:** `python3 demo_queries.py`. Current volumes: ~3.58M titers; ~102k sequences
  (~38%); ~3.85k locations (~99.98% resolved).
- **Git operations need the sandbox disabled** (writes to `.git` are blocked).

## Pipeline

Each stage writes a CSV that `load_duckdb.py` loads.

| Stage | Does | Cost |
|---|---|---|
| `build_db.py` | hidb5 → antigen/serum/titer_table/titer | pure python, ~5 s |
| `build_locations.py` | locationdb → location | pure python, ~1 s |
| `build_sequences.py` | seqdb → match (ae_backend, parallel, **incremental cache** keyed on a stable natural key) | ~4 min cold / ~seconds warm |
| `build_clades.py` | clades.json → clade (ae_backend) | ~4 s |
| `load_duckdb.py` | CSVs → `serology.duckdb` + Parquet; tables antigen/serum/titer_table/titer/sequence/location; views `titer_flat`, `antigen_sequence` | ~7 s |
| `refresh.sh` | orchestrates via two hashes: **DATA** (hidb5+seqdb+locationdb) and **CLADE** (clades.json); runs only what changed | — |

## Graduation objectives

### Data fidelity (highest scientific value)

1. **Reconcile against source `.ace` charts.** The DB is built from hidb5's *merged*
   view, and ~0.9% of titer cells (~32.7k) are dropped where a table has more titer rows
   than registered hidb indices. Decide: parse `.ace` directly (whocc-tables) for
   canonical per-table provenance, or reconcile just the clipped cells against `.ace`.
   Report exact recovered counts.
2. **Validate scientific conventions** against authoritative acmacs/Racmacs computations
   on a sample: log-titer = log2(titer/10), GMT, and threshold-censoring. Censored titers
   are stored with a `kind` flag but downstream stats treat them simply — decide the
   correct (censored-aware) handling.
3. **Sequence coverage (~38%)** is mostly genuine (unsequenced antigens); the ae matcher
   agrees with a naive join. The lab_id fallback is **inert** against the current seqdb
   (empty lab_id lists — already verified; do **not** chase it). A real lever is exposing
   ae's name-indexed `select_by_name` / `populate_from_seqdb` to Python (needs a pybind
   binding + ae rebuild) — evaluate whether it's worth it.

### Engineering hardening

4. **Test suite (none exists).** Add data-quality assertions (row counts, join integrity,
   no orphan FKs, titer/log-titer ranges, resolution rates) + regression tests on a small
   fixture. Do this **early** so later changes are safe.
5. **Schema constraints + stable keys.** Add enforced PK/FK where DuckDB allows; replace
   positional IDs (e.g. `h3:1042`, unstable across regenerations) with content-based
   natural keys — prerequisite for true incremental and the eventual Postgres store.
6. **Reproducibility.** Pin deps (`pyproject`/requirements), document the exact env,
   consider a smoke-test CI. The toolchain is currently fragile.
7. **(Only if the user explicitly wants it now) the service.** The README describes a
   Postgres upload/query/download service — a larger build. Otherwise keep the schema
   Postgres-portable.

## Hard constraints

- **The repo is PUBLIC and must never contain WHO data** — no strain names, serum/lab
  IDs, GISAID accessions, sequence strings, or titer values in tracked files or commit
  messages. A pre-commit hook enforces this; activate it once with
  `git config core.hooksPath hooks`. All data stays in the gitignored `$SERO_OUT`.
- **Push only to `drserajames`, never `acorg`.** Commit/push only when the user asks.
- **Keep it portable** — no hardcoded `/Users/...` paths; defaults resolve relative to the
  repo (siblings `acmacs-data/`, `ae/`), overridable via `ACMACS_DATA` / `AE_BUILD` /
  `SERO_OUT`.
- **`whocc-tables/` (the `.ace` source) and `fludata/` are read-only** — don't modify.

## How to work

- **Verify empirically at each step:** run the pipeline, check row counts, run
  `demo_queries.py` and your tests. Report fidelity numbers — don't just assert.
- **Keep the prototype working** after every change (it's used for analysis now).
- **Before large/architectural changes** (parse `.ace` vs stay hidb-derived; content-key
  migration; starting the service), write a short plan and get user approval — several of
  these are the user's call.
- **Update the README** as behavior changes.

## First steps

1. Read `README.md` end to end.
2. Reproduce: run `./refresh.sh` and `demo_queries.py`; confirm `ae_backend` imports and
   the numbers match the volumes above.
3. Propose a prioritized graduation plan for approval. Suggested order: test suite (#4)
   first, then the `.ace` fidelity reconciliation (#1) and scientific validation (#2) —
   those remove the biggest caveats.
