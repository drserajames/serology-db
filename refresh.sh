#!/usr/bin/env bash
# refresh.sh — rebuild the serology-db prototype, running only the stale stages.
#
#   ./refresh.sh           rebuild whatever changed (data and/or clades)
#   ./refresh.sh --force   rebuild everything regardless of hashes
#   ./refresh.sh --no-seq  titers only; skip the (slow, ~4 min) sequence match
#   ./refresh.sh --pull    first pull fresh hidb5 from upstream (needs WHO CC
#                          SSH access via acmacs-data/hidb5-update-download)
#
# Two independent cadences, hashed separately so each re-runs only when needed:
#   DATA  = hidb5 + seqdb + locationdb   -> titers (build_db) + locations
#           (build_locations, ~1 s) + match (build_sequences).  The match is
#           INCREMENTAL: ~4 min on first build, then ~5-30 s (only new /
#           newly-matchable keys).  Changes ~1-2 weeks.
#   CLADE = clades.json                          -> clades (build_clades, ~4 s).
#           A clades.json edit re-derives clades WITHOUT re-running the match.
# load_duckdb re-assembles the DB from CSVs whenever any stage ran.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ACMACS_DATA="${ACMACS_DATA:-/Users/sarahjames/AC/eu/acmacs-data}"
# Output lives outside this code dir (large, WHO-derived); gitignored in
# acmacs-data. Override with SERO_OUT. The python scripts read the same default.
OUT="${SERO_OUT:-$ACMACS_DATA/serology-db}"
export ACMACS_DATA SERO_OUT="$OUT"
MANIFEST="$OUT/.manifest"
DATA_FILES=("$ACMACS_DATA"/hidb5.{h1,h3,b}.json.xz
            "$ACMACS_DATA"/seqdb-{h1,h3,b}.v4.json.xz
            "$ACMACS_DATA"/locationdb.json.xz)
CLADES_FILE="$ACMACS_DATA/clades.json"

FORCE=0; PULL=0; SEQ=1
for arg in "$@"; do
  case "$arg" in
    --force)  FORCE=1 ;;
    --pull)   PULL=1 ;;
    --no-seq) SEQ=0 ;;
    -h|--help) sed -n 's/^# //p' "$0"; exit 0 ;;
    *) echo "unknown option: $arg (try --help)" >&2; exit 2 ;;
  esac
done

if [[ $PULL == 1 ]]; then
  if [[ -x "$ACMACS_DATA/hidb5-update-download" ]]; then
    echo ">> pulling fresh hidb5 from upstream (needs WHO CC SSH access)…"
    "$ACMACS_DATA/hidb5-update-download"
  else
    echo "!! $ACMACS_DATA/hidb5-update-download not found/executable; skipping pull" >&2
  fi
fi

for f in "${DATA_FILES[@]}"; do
  [[ -f "$f" ]] || { echo "!! missing data source: $f" >&2; exit 1; }
done
mkdir -p "$OUT"

cur_data="$(shasum -a 256 "${DATA_FILES[@]}" | shasum -a 256 | awk '{print $1}')"
cur_clade="$([[ -f "$CLADES_FILE" ]] && shasum -a 256 "$CLADES_FILE" | awk '{print $1}' || echo none)"
prev_data="$(awk '/^DATA/{print $2}'  "$MANIFEST" 2>/dev/null || true)"
prev_clade="$(awk '/^CLADE/{print $2}' "$MANIFEST" 2>/dev/null || true)"
prev_seq="$(awk '/^SEQ/{print $2}'    "$MANIFEST" 2>/dev/null || echo 0)"

# Decide which stages are stale.
titer_stale=0; seq_stale=0; clade_stale=0
[[ $FORCE == 1 || "$cur_data" != "$prev_data" || ! -f "$OUT/serology.duckdb" ]] && titer_stale=1
if [[ $SEQ == 1 ]]; then
  [[ $FORCE == 1 || "$cur_data" != "$prev_data" || "$prev_seq" != 1 ]] && seq_stale=1
  [[ $FORCE == 1 || "$cur_clade" != "$prev_clade" || $seq_stale == 1 ]] && clade_stale=1
fi

if [[ $titer_stale == 0 && $seq_stale == 0 && $clade_stale == 0 ]]; then
  echo "== up to date — nothing changed (data $cur_data / clade $cur_clade)."
  exit 0
fi

cd "$HERE"
if [[ $titer_stale == 1 ]]; then
  echo ">> titers (hidb5)…";           python3 build_db.py
  echo ">> locations (locationdb)…";   python3 build_locations.py || echo ">> locations skipped" >&2
fi
if [[ $seq_stale == 1 ]]; then
  echo ">> sequence match (incremental; ~4 min first build, else ~5-30 s)…"
  python3 build_sequences.py || echo ">> match skipped (ae unavailable); titer DB intact" >&2
fi
if [[ $clade_stale == 1 ]]; then
  echo ">> clades (clades.json, ~4 s)…"
  python3 build_clades.py || echo ">> clades skipped (ae unavailable)" >&2
fi
python3 load_duckdb.py

# Record what the rebuilt artifacts actually reflect.
new_seq=0;   [[ $SEQ == 1 && -f "$OUT/csv/match.csv" ]] && new_seq=1
new_clade="$prev_clade"
[[ $clade_stale == 1 && -f "$OUT/csv/clade.csv" ]] && new_clade="$cur_clade"
{ echo "DATA $cur_data"; echo "CLADE $new_clade"; echo "SEQ $new_seq"; } > "$MANIFEST"
echo "== rebuilt (data=$cur_data clade=$new_clade seq=$new_seq)."
