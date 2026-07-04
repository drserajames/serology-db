#!/usr/bin/env python3
"""Recover CLIPPED titer cells from the canonical source .ace charts (optional stage).

build_db drops titer cells that fall outside a hidb5 table's registered antigen/
serum indices — antigens/sera hidb5 itself could not resolve to its identity DB
(~0.9% of cells; ae's own hidb reader drops them identically). Those readings are
real, though: the raw source chart in whocc-tables/ carries the full pre-resolution
antigen/serum list, and its matrix aligns 1:1 with hidb's (verified below).

This stage maps each clipped table to its source .ace, and — ONLY when the source's
registered cells are byte-identical to hidb's matrix (strict alignment, so we never
attach a wrong strain identity) — emits the previously-dropped rows/cols as new
`antigen`/`serum` rows plus their `titer` cells, tagged provenance `ace_recovered`.

  hidb table (rows 0..len(a)-1 registered, len(a)..N-1 clipped)
     -> source .ace (antigens 0..N-1)   [require registered cells == .ace cells]
     -> recovered antigens = .ace[len(a):] , recovered sera = .ace_sera[len(s):]
     -> recovered titers   = cells touching a clipped row or column

Recovered antigens/sera get CONTENT-DERIVED ids (`{sub}:ra:{hash}` / `{sub}:rs:{hash}`
of the natural key name+passage+reassortant), so the same unresolved strain tested
across several tables collapses to one row and ids never collide with hidb's
positional `{sub}:{i}`. Recovered titers keep the ORIGINAL `tab_id`, so they attach
to the existing titer_table row. (Recovered antigens are deduped among themselves
but NOT merged into hidb antigens — these are exactly the strains hidb could not
identity-resolve, so a name-based merge would risk the very errors hidb avoided.)

Input : acmacs-data/hidb5.{h1,h3,b}.json.xz , whocc-tables/*.ace
Output: out/csv/recovered_{antigen,serum,titer}.csv
Env (auto-defaulted): ACMACS_DATA, WHOCC_TABLES, PYTHONPATH->AE_BUILD.
OPTIONAL: if ae_backend can't load or whocc-tables is absent, this is skipped and
the titer DB is unaffected (recovered_*.csv simply aren't produced).
"""
import csv
import glob
import hashlib
import json
import lzma
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ACMACS_DATA = os.environ.get(
    "ACMACS_DATA", os.path.normpath(os.path.join(HERE, os.pardir, "acmacs-data")))
WHOCC = os.environ.get(
    "WHOCC_TABLES", os.path.normpath(os.path.join(HERE, os.pardir, "whocc-tables")))
AE_BUILD = os.environ.get(
    "AE_BUILD", os.path.normpath(os.path.join(HERE, os.pardir, "ae", "build")))
OUT = os.environ.get("SERO_OUT", os.path.join(ACMACS_DATA, "serology-db"))
CSV_DIR = os.environ.get("SERO_CSV_DIR", os.path.join(OUT, "csv"))
sys.path.insert(0, AE_BUILD)

SUBTYPES = ["h1", "h3", "b"]
SUBTYPE_TAG = {"h1": "A(H1N1)", "h3": "A(H3N2)", "b": "B"}
_NAME_PARSE = None   # set to ae_backend.virus.name_parse in main() (post-import)
SUBFAM = {"h1": ["h1pdm", "h1seas"], "h3": ["h3"], "b": ["b", "bvic", "byam"]}
ASSAY = {"HI": "hi", "FRA": "fra", "FOCUS REDUCTION": "fr", "MN": "mn",
         "HINT": "hint", "PRN": "prn", "PRNT": "prn"}
_titer_re = re.compile(r"^([<>~]?)(\d+)")


def load_json_xz(path):
    return json.loads(lzma.open(path, "rb").read().decode("utf-8").replace("\\U", "\\u"))


def parse_titer(raw):
    """(kind, value, log_titer); '*'/missing -> None. Mirrors build_db.parse_titer."""
    if raw is None:
        return None
    s = str(raw).strip()
    if s in ("", "*"):
        return None
    m = _titer_re.match(s)
    if not m:
        return ("other", None, None)
    import math
    pfx, num = m.group(1), int(m.group(2))
    kind = {"<": "lt", ">": "gt", "~": "num", "": "num"}[pfx]
    log_titer = math.log2(num / 10.0) if num > 0 else None
    return (kind, num, log_titer)


def _attr(obj, name):
    """Read an ae accessor that may be a method or a property; '' on failure."""
    try:
        v = getattr(obj, name)
        v = v() if callable(v) else v
        return "" if v is None else str(v)
    except Exception:
        return ""


def candidate_charts(sub, tb):
    assay = ASSAY.get(str(tb.get("A", "")).upper())
    rbc = str(tb.get("r", "") or "").lower().replace(" ", "-")
    lab = str(tb.get("l", "") or "").lower()
    date = str(tb.get("D", ""))
    if not (assay and lab and date):
        return []
    hits = []
    for fam in SUBFAM.get(sub, []):
        for d in glob.glob(os.path.join(WHOCC, f"{fam}-{assay}-{rbc}-{lab}")):
            hits += glob.glob(os.path.join(d, f"*-{date}.ace"))
            hits += glob.glob(os.path.join(d, f"*-{date}_*.ace"))
    return sorted(set(hits))


def aligned_chart(c3, sub, tb):
    """Return the source Chart iff dims match AND every registered cell is identical."""
    na_reg, ns_reg = len(tb["a"]), len(tb["s"])
    ncol = max((len(r) for r in tb["t"]), default=0)
    for path in candidate_charts(sub, tb):
        ch = c3.Chart(path)
        # full matrix must align: antigens == rows, sera == widest row
        if ch.number_of_antigens() != len(tb["t"]) or ch.number_of_sera() != ncol:
            continue
        titers = ch.titers()
        ok = True
        for r in range(na_reg):
            row = tb["t"][r]
            for c in range(ns_reg):
                if str(row[c]) != str(titers.titer(r, c)):
                    ok = False
                    break
            if not ok:
                break
        if ok:
            return ch
    return None


def _rec_id(sub, prefix, *parts):
    """Content-derived id from a natural key — stable across tables and rebuilds,
    so the same unresolved strain (tested in several tables) collapses to one row.
    Namespace `{sub}:{prefix}:{hash}` can't collide with hidb's `{sub}:{int}`."""
    h = hashlib.md5("|".join([sub, *(str(p) for p in parts)]).encode()).hexdigest()
    return f"{sub}:{prefix}:{h[:12]}"


def _loc_year(name):
    """Location + year via ae's canonical virus-name parser (matches locationdb).
    The .ace Antigen has no location/year accessor, so parse the full name. Returns
    ('','') for the ~17% of reference/control names the parser can't split."""
    if not (_NAME_PARSE and name):
        return "", ""
    try:
        parts = _NAME_PARSE(name).parts
        return str(parts.location), str(parts.year)
    except Exception:
        return "", ""


def antigen_row(sub, ag):
    """Return (ag_id, row-fields) for a recovered .ace antigen."""
    name, passage, reass = _attr(ag, "name"), _attr(ag, "passage"), _attr(ag, "reassortant")
    location, year = _loc_year(name)
    ag_id = _rec_id(sub, "ra", name, passage, reass)
    return ag_id, [ag_id, sub, name, SUBTYPE_TAG.get(sub, ""),
                   "", location, "", year, passage, reass, _attr(ag, "date")]


def serum_row(sub, sr):
    """Return (sr_id, row-fields) for a recovered .ace serum."""
    sid, name = _attr(sr, "serum_id"), _attr(sr, "name")
    passage, reass = _attr(sr, "passage"), _attr(sr, "reassortant")
    sr_id = _rec_id(sub, "rs", sid, name, passage, reass)
    return sr_id, [sr_id, sub, sid, name, SUBTYPE_TAG.get(sub, ""),
                   _attr(sr, "lineage"), "", "", "", passage, _attr(sr, "serum_species")]


def main():
    try:
        import ae_backend
    except Exception as err:
        print(f"!! ae_backend unavailable ({err}); skipping clipped recovery "
              f"(titer DB unaffected).", file=sys.stderr)
        return 1
    if not os.path.isdir(WHOCC):
        print(f"!! whocc-tables not found at {WHOCC}; skipping clipped recovery.",
              file=sys.stderr)
        return 1
    c3 = ae_backend.chart_v3
    global _NAME_PARSE
    _NAME_PARSE = ae_backend.virus.name_parse
    os.makedirs(CSV_DIR, exist_ok=True)

    ag_f = open(os.path.join(CSV_DIR, "recovered_antigen.csv"), "w", newline="")
    sr_f = open(os.path.join(CSV_DIR, "recovered_serum.csv"), "w", newline="")
    ti_f = open(os.path.join(CSV_DIR, "recovered_titer.csv"), "w", newline="")
    ag_w, sr_w, ti_w = csv.writer(ag_f), csv.writer(sr_f), csv.writer(ti_f)
    ag_w.writerow(["ag_id", "subtype", "name", "virus_type", "lineage",
                   "location", "isolation", "year", "passage", "reassortant",
                   "collection_date"])
    sr_w.writerow(["sr_id", "subtype", "serum_id", "name", "virus_type",
                   "lineage", "location", "isolation", "year", "passage", "species"])
    ti_w.writerow(["tab_id", "ag_id", "sr_id", "titer_raw", "titer_kind",
                   "titer_value", "log_titer"])

    seen_ag, seen_sr = set(), set()
    n_tab = n_rec_ag = n_rec_sr = n_rec_ti = 0
    n_clip_tables = n_aligned = 0
    for sub in SUBTYPES:
        path = os.path.join(ACMACS_DATA, f"hidb5.{sub}.json.xz")
        if not os.path.exists(path):
            continue
        for ti, tb in enumerate(load_json_xz(path)["t"]):
            ncol = max((len(r) for r in tb["t"]), default=0)
            if not (len(tb["t"]) > len(tb["a"]) or ncol > len(tb["s"])):
                continue
            n_clip_tables += 1
            ch = aligned_chart(c3, sub, tb)
            if ch is None:
                continue
            n_aligned += 1
            tab_id = f"{sub}:{ti}"
            na_reg, ns_reg = len(tb["a"]), len(tb["s"])

            # map matrix index -> ag_id / sr_id (registered use hidb id; clipped new)
            def ag_of(r):
                if r < na_reg:
                    return f"{sub}:{tb['a'][r]}"
                rid, fields = antigen_row(sub, ch.antigen(r))
                if rid not in seen_ag:
                    seen_ag.add(rid)
                    ag_w.writerow(fields)
                return rid

            def sr_of(c):
                if c < ns_reg:
                    return f"{sub}:{tb['s'][c]}"
                rid, fields = serum_row(sub, ch.serum(c))
                if rid not in seen_sr:
                    seen_sr.add(rid)
                    sr_w.writerow(fields)
                return rid

            before_ag, before_sr = len(seen_ag), len(seen_sr)
            for r, row in enumerate(tb["t"]):
                for c, cell in enumerate(row):
                    if r < na_reg and c < ns_reg:
                        continue  # registered cell — already in build_db
                    parsed = parse_titer(cell)
                    if parsed is None:
                        continue
                    kind, val, lg = parsed
                    ti_w.writerow([tab_id, ag_of(r), sr_of(c), cell, kind, val, lg])
                    n_rec_ti += 1
            n_rec_ag += len(seen_ag) - before_ag
            n_rec_sr += len(seen_sr) - before_sr
            n_tab += 1

    for f in (ag_f, sr_f, ti_f):
        f.close()
    print(f"  clipped tables={n_clip_tables} strict-aligned={n_aligned} "
          f"recovered from {n_tab} tables")
    print(f"TOTAL  recovered antigens={n_rec_ag} sera={n_rec_sr} titers={n_rec_ti} "
          f"->  out/csv/recovered_*.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
