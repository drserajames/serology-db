#!/usr/bin/env python3
"""Build a local DuckDB + Parquet prototype from the hidb5 titer databases.

Source: acmacs-data/hidb5.{h1,h3,b}.json.xz  (already-normalised antigen/serum
identity + the raw WHO CC titer matrices).  We reshape the per-table titer
matrices into a tidy long `titer` fact plus antigen / serum / table dimensions
as CSV.  (HA sequences + clades are attached separately by build_sequences.py,
which reuses ae_backend's canonical seqdb matcher; load_duckdb.py then loads it
all and exports Parquet.)

  hidb5 (nested, per-subtype) ->  CSV  ->  [DuckDB (typed)  ->  Parquet]

NOTE: contains real WHO CC serology data. Keep local; do not commit/push.
"""
import csv, json, lzma, math, os, re, sys

ACMACS_DATA = os.environ.get("ACMACS_DATA", "/Users/sarahjames/AC/eu/acmacs-data")
# Output lives OUTSIDE this code dir (large, WHO-derived) — defaults into
# acmacs-data (gitignored there). Override with SERO_OUT.
OUT = os.environ.get("SERO_OUT", os.path.join(ACMACS_DATA, "serology-db"))
CSV_DIR = os.environ.get("SERO_CSV_DIR", os.path.join(OUT, "csv"))
SUBTYPES = ["h1", "h3", "b"]

_titer_re = re.compile(r"^([<>~]?)(\d+)")


def load_json_xz(path):
    """Load an acmacs .json.xz (hidb5 or seqdb), normalising the non-standard
    \\U#### CJK escapes (e.g. Chinese location names) to standard JSON \\u####."""
    raw = lzma.open(path, "rb").read().decode("utf-8")
    raw = raw.replace("\\U", "\\u")
    return json.loads(raw)


def parse_titer(raw):
    """Return (kind, value, log_titer). kind in num/lt/gt/other; '*' -> None."""
    if raw is None:
        return None
    s = str(raw).strip()
    if s in ("", "*"):
        return None
    m = _titer_re.match(s)
    if not m:
        return ("other", None, None)
    pfx, num = m.group(1), int(m.group(2))
    kind = {"<": "lt", ">": "gt", "~": "num", "": "num"}[pfx]
    log_titer = math.log2(num / 10.0) if num > 0 else None
    return (kind, num, log_titer)


def name_from(rec):
    """Reconstruct a strain name O/i/y, tolerating heterogeneous records."""
    iso = rec.get("i", "")
    if "/" in str(iso):           # already a full name (older serum records)
        return iso
    loc, yr = rec.get("O", ""), rec.get("y", "")
    return "/".join(p for p in (loc, iso, yr) if p)


def first(x):
    return x[0] if isinstance(x, list) and x else (x or None)


def main():
    os.makedirs(CSV_DIR, exist_ok=True)

    ag_f = open(os.path.join(CSV_DIR, "antigen.csv"), "w", newline="")
    sr_f = open(os.path.join(CSV_DIR, "serum.csv"), "w", newline="")
    tb_f = open(os.path.join(CSV_DIR, "titer_table.csv"), "w", newline="")
    ti_f = open(os.path.join(CSV_DIR, "titer.csv"), "w", newline="")
    ag_w, sr_w, tb_w, ti_w = (csv.writer(f) for f in (ag_f, sr_f, tb_f, ti_f))

    ag_w.writerow(["ag_id", "subtype", "name", "virus_type", "lineage",
                   "location", "isolation", "year", "passage", "reassortant",
                   "collection_date"])
    sr_w.writerow(["sr_id", "subtype", "serum_id", "name", "virus_type",
                   "lineage", "location", "isolation", "year", "passage", "species"])
    tb_w.writerow(["tab_id", "subtype", "assay", "rbc", "lab", "virus",
                   "table_date"])
    ti_w.writerow(["tab_id", "ag_id", "sr_id", "titer_raw", "titer_kind",
                   "titer_value", "log_titer"])

    n_ag = n_sr = n_tb = n_ti = n_clip = 0
    for sub in SUBTYPES:
        path = os.path.join(ACMACS_DATA, f"hidb5.{sub}.json.xz")
        if not os.path.exists(path):
            print(f"  skip {sub}: {path} missing", file=sys.stderr); continue
        d = load_json_xz(path)
        A, S, T = d["a"], d["s"], d["t"]

        for i, a in enumerate(A):
            ag_w.writerow([f"{sub}:{i}", sub, name_from(a), a.get("V", ""),
                           a.get("L", ""), a.get("O", ""), a.get("i", ""),
                           a.get("y", ""), a.get("P", ""), a.get("R", ""),
                           first(a.get("D"))])
        for i, s in enumerate(S):
            sr_w.writerow([f"{sub}:{i}", sub, s.get("I", ""), name_from(s),
                           s.get("V", ""), s.get("L", ""), s.get("O", ""),
                           s.get("i", ""), s.get("y", ""), s.get("P", ""),
                           s.get("s", "")])
        n_ag += len(A); n_sr += len(S)

        for ti, tb in enumerate(T):
            tab_id = f"{sub}:{ti}"
            raw_d = str(tb.get("D", ""))
            tdate = (f"{raw_d[0:4]}-{raw_d[4:6]}-{raw_d[6:8]}"
                     if len(raw_d) == 8 and raw_d.isdigit() else None)
            tb_w.writerow([tab_id, sub, tb.get("A", ""), tb.get("r", ""),
                           tb.get("l", ""), tb.get("V", ""), tdate])
            ag_idx, sr_idx, matrix = tb["a"], tb["s"], tb["t"]
            # A few tables carry more titer rows/cols than registered hidb
            # indices (antigens/sera not in the identity DB). Clip to the
            # indexed extent and tally what was dropped.
            for r, row in enumerate(matrix):
                if r >= len(ag_idx):
                    n_clip += sum(1 for c in row if parse_titer(c))
                    continue
                ag_id = f"{sub}:{ag_idx[r]}"
                for c, cell in enumerate(row):
                    if c >= len(sr_idx):
                        if parse_titer(cell):
                            n_clip += 1
                        continue
                    parsed = parse_titer(cell)
                    if parsed is None:        # missing '*' -> omit
                        continue
                    kind, val, lg = parsed
                    ti_w.writerow([tab_id, ag_id, f"{sub}:{sr_idx[c]}",
                                   cell, kind, val, lg])
                    n_ti += 1
        n_tb += len(T)
        print(f"  {sub}: {len(A)} antigens, {len(S)} sera, {len(T)} tables")

    for f in (ag_f, sr_f, tb_f, ti_f):
        f.close()
    print(f"TOTAL  antigens={n_ag}  sera={n_sr}  tables={n_tb}  "
          f"titers={n_ti}  clipped_cells={n_clip}")
    return n_ti


if __name__ == "__main__":
    main()
