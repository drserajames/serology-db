#!/usr/bin/env python3
"""Derive clades for matched sequences — the fast, clades.json-dependent stage.

find_clades() over the ENTIRE seqdb takes ~0.16 s, so re-deriving clades when
clades.json changes is near-instant and does NOT require re-running the (slow)
name match.  Clades depend only on the sequence, so we key by seq_id.

Input : out/csv/match.csv           (ag_id -> seq_id, from build_sequences.py)
Output: out/csv/clade.csv           seq_id -> clade, clade_path

Env (auto-defaulted): PYTHONPATH→ae/build, SEQDB_V4, LOCDB_V2, AC_CLADES_JSON_V2.
"""
import csv, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ACMACS_DATA = os.environ.get("ACMACS_DATA", "/Users/sarahjames/AC/eu/acmacs-data")
AE_BUILD = os.environ.get("AE_BUILD", "/Users/sarahjames/AC/eu/ae/build")
OUT = os.environ.get("SERO_OUT", os.path.join(ACMACS_DATA, "serology-db"))
CSV_DIR = os.environ.get("SERO_CSV_DIR", os.path.join(OUT, "csv"))
os.environ.setdefault("SEQDB_V4", ACMACS_DATA)
os.environ.setdefault("LOCDB_V2", os.path.join(ACMACS_DATA, "locationdb.json.xz"))
os.environ.setdefault("AC_CLADES_JSON_V2", os.path.join(ACMACS_DATA, "clades.json"))
sys.path.insert(0, AE_BUILD)

SUBTYPE_TAG = {"h1": "A(H1N1)", "h3": "A(H3N2)", "b": "B"}
CLADES_FILE = os.environ["AC_CLADES_JSON_V2"]
MATCH = os.path.join(CSV_DIR, "match.csv")


def main():
    if not os.path.exists(MATCH):
        print(f"!! {MATCH} missing; run build_sequences.py first.", file=sys.stderr)
        return 1
    try:
        import ae_backend
    except Exception as err:
        print(f"!! ae_backend unavailable ({err}); skipping clades.", file=sys.stderr)
        return 1

    # which seq_ids do we need clades for, per subtype
    need = {}
    with open(MATCH) as f:
        for r in csv.DictReader(f):
            need.setdefault(r["subtype"], set()).add(r["seq_id"])

    out = open(os.path.join(CSV_DIR, "clade.csv"), "w", newline="")
    w = csv.writer(out)
    w.writerow(["seq_id", "clade", "clade_path"])

    n = 0
    for sub, seq_ids in need.items():
        tag = SUBTYPE_TAG.get(sub)
        if tag is None:
            continue
        seqdb = ae_backend.seqdb.for_subtype(tag)
        sel = seqdb.select_all()
        sel.find_clades(CLADES_FILE)        # once over the whole subtype (~0.2 s)
        for i in range(len(sel)):
            ref = sel[i]
            sid = ref.seq_id()
            if sid not in seq_ids:
                continue
            clades = list(ref.clades or [])
            w.writerow([sid, clades[-1] if clades else "", " ".join(clades)])
            n += 1
        print(f"  {sub}: clades for {len(seq_ids)} matched sequences")
    out.close()
    print(f"TOTAL  {n} sequence clades  ->  out/csv/clade.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
