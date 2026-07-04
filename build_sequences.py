#!/usr/bin/env python3
"""Match antigens to HA sequences via ae_backend's canonical matcher — incrementally.

Reuses ae's seqdb matcher (as ae/py/ae/report/geographic.py does).  filter_name()
is CPU-bound, does NOT release the GIL, and names are ~94% distinct, so we fan the
distinct keys across processes.  A full match is ~4 min — but most keys don't
change between the biweekly refreshes, so we cache results by the STABLE natural
key (subtype, name, reassortant, passage) and only (re)match what's actually new
or newly matchable:

  * cached HIT  -> reused (sequences aren't removed, so a name hit stays valid)
  * new key     -> always matched
  * cached MISS -> re-checked only when seqdb changed, and (by default) only for
                   recent strains — old unmatched strains are ~never newly
                   sequenced, so re-trying 132k of them every fortnight is waste.

Typical incremental cost: ~5 s (seqdb unchanged) to ~30 s (seqdb grew), vs ~240 s.

Input : out/csv/antigen.csv
State : out/csv/match_cache.csv    persistent per-key cache (hits AND misses)
        out/csv/.match_meta        seqdb hash + match mode of the last run
Output: out/csv/match.csv          ag_id -> seq_id + aa (derived from the cache)

Flags: --with-passage      passage-specific HA (own cache namespace)
       --index             opt-in fast path: build a name index once and dict-look
                           up each key instead of a full seqdb scan per key. Cuts a
                           cold match from ~4 min to ~seconds; byte-identical result
                           (falls back to filter_name for the ~2-3% ambiguous names).
       --rematch-all       ignore the cache, match every current key
       --rematch-misses    on seqdb change, re-check ALL misses (not just recent)
       --recheck-years N   recency window for miss re-checks (default 3)
Env (auto-defaulted): PYTHONPATH→ae/build, SEQDB_V4, LOCDB_V2.  Optional: if
ae_backend can't load, the titer DB is unaffected.
"""
import csv, hashlib, os, sys
from multiprocessing import Pool, cpu_count

HERE = os.path.dirname(os.path.abspath(__file__))
# Defaults assume acmacs-data and ae are siblings of this repo; override via env.
ACMACS_DATA = os.environ.get("ACMACS_DATA", os.path.normpath(os.path.join(HERE, os.pardir, "acmacs-data")))
AE_BUILD = os.environ.get("AE_BUILD", os.path.normpath(os.path.join(HERE, os.pardir, "ae", "build")))
OUT = os.environ.get("SERO_OUT", os.path.join(ACMACS_DATA, "serology-db"))
CSV_DIR = os.environ.get("SERO_CSV_DIR", os.path.join(OUT, "csv"))
CACHE = os.path.join(CSV_DIR, "match_cache.csv")
META = os.path.join(CSV_DIR, ".match_meta")
os.environ.setdefault("SEQDB_V4", ACMACS_DATA)
os.environ.setdefault("LOCDB_V2", os.path.join(ACMACS_DATA, "locationdb.json.xz"))

SUBTYPE_TAG = {"h1": "A(H1N1)", "h3": "A(H3N2)", "b": "B"}
WITH_PASSAGE = "--with-passage" in sys.argv
USE_INDEX = "--index" in sys.argv          # opt-in pure-Python name-index fast path
REMATCH_ALL = "--rematch-all" in sys.argv
REMATCH_MISSES = "--rematch-misses" in sys.argv
RECHECK_YEARS = 3
if "--recheck-years" in sys.argv:
    RECHECK_YEARS = int(sys.argv[sys.argv.index("--recheck-years") + 1])
NPROC = min(cpu_count(), 10)
csv.field_size_limit(1 << 24)

_SDB = None


def _init(tag):
    global _SDB
    sys.path.insert(0, AE_BUILD)
    import ae_backend
    _SDB = ae_backend.seqdb.for_subtype(tag)


def _work(keys):
    """keys: (name, reassortant, passage) -> {key: (seq_id, aa, nuc_len)} (hits only)."""
    out = {}
    for name, reassortant, passage in keys:
        sel = _SDB.select_all().filter_name(name=name, reassortant=reassortant,
                                            passage=passage)
        if len(sel):
            ref = sel[0]
            out[(name, reassortant, passage)] = (ref.seq_id(), str(ref.aa),
                                                 len(str(ref.nuc)))
    return out


def _match_indexed(tag, triples):
    """Opt-in (--index) fast path, byte-identical to _work.

    The slow path's cost is rebuilding the whole `select_all()` selection for every
    key (~2 ms each, memory-bandwidth-bound). Instead, iterate the seqdb ONCE into a
    name -> [ref] map and dict-lookup each key:
      * name absent      -> miss (instant)
      * one seq for name -> that seq (filter_name would return the same single
                            survivor; no ranking happens for size==1)
      * >1 seq for name  -> fall back to select_all().filter_name for exact ranking
                            (only ~2-3% of keys; keeps the result identical)
    So misses + single-candidate names (97-98% of keys) skip select_all() entirely.
    Runs single-process (the lookups are trivial once the index is built)."""
    import ae_backend
    sdb = ae_backend.seqdb.for_subtype(tag)
    s = sdb.select_all()
    by_name = {}
    for i in range(len(s)):
        by_name.setdefault(s[i].name(), []).append(i)   # .name() == entry name (no subtype prefix)
    prefix = f"{tag}/"
    out = {}
    for (name, reassortant, passage) in triples:
        nm = name[len(prefix):] if name.startswith(prefix) else name
        cands = by_name.get(nm)
        if not cands:
            continue
        if len(cands) == 1:
            ref = s[cands[0]]
        else:
            sel = sdb.select_all().filter_name(name=name, reassortant=reassortant,
                                               passage=passage)
            if not len(sel):
                continue
            ref = sel[0]
        out[(name, reassortant, passage)] = (ref.seq_id(), str(ref.aa),
                                             len(str(ref.nuc)))
    return out


def chunked(seq, n):
    k, m = divmod(len(seq), n)
    return [seq[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n)]


def seqdb_hash():
    h = hashlib.sha256()
    for sub in ("h1", "h3", "b"):
        p = os.path.join(ACMACS_DATA, f"seqdb-{sub}.v4.json.xz")
        if os.path.exists(p):
            with open(p, "rb") as f:
                for blk in iter(lambda: f.read(1 << 20), b""):
                    h.update(blk)
    return h.hexdigest()


def load_cache():
    cache = {}                     # (subtype,name,reass,passage) -> (seq_id,aal,nucl,aa)
    if os.path.exists(CACHE):
        with open(CACHE) as f:
            for r in csv.DictReader(f):
                cache[(r["subtype"], r["name"], r["reassortant"], r["passage"])] = (
                    r["seq_id"], r["aa_length"], r["nuc_length"], r["aa"])
    meta = {}
    if os.path.exists(META):
        with open(META) as f:
            for line in f:
                if " " in line:
                    k, v = line.rstrip("\n").split(" ", 1)
                    meta[k] = v
    return cache, meta


def main():
    try:
        sys.path.insert(0, AE_BUILD)
        import ae_backend  # noqa: F401
    except Exception as err:
        print(f"!! ae_backend unavailable ({err}); skipping sequences "
              f"(titer DB unaffected).", file=sys.stderr)
        return 1

    # current antigens -> (ag_id, key, year), grouped by subtype
    rows = []
    with open(os.path.join(CSV_DIR, "antigen.csv")) as f:
        for r in csv.DictReader(f):
            passage = r["passage"] if WITH_PASSAGE else ""
            key = (r["subtype"], r["name"], r["reassortant"], passage)
            y = r["year"] or (r["collection_date"][:4] if r["collection_date"] else "")
            rows.append((r["ag_id"], key, y))
    cur_keys = {k for _, k, _ in rows}
    year_of = {}
    for _, k, y in rows:
        year_of.setdefault(k, y)

    cache, meta = load_cache()
    cur_hash = seqdb_hash()
    mode = "passage" if WITH_PASSAGE else "name"
    mode_changed = meta.get("mode") != mode
    seqdb_changed = meta.get("seqdb") != cur_hash

    def as_int(y):
        try:
            return int(y)
        except (TypeError, ValueError):
            return -1
    max_year = max((as_int(y) for y in year_of.values()), default=0)
    threshold = max_year - RECHECK_YEARS

    if REMATCH_ALL or mode_changed or not cache:
        to_match = set(cur_keys)
        reused = 0
    else:
        new = {k for k in cur_keys if k not in cache}
        recheck = set()
        if seqdb_changed:
            for k in cur_keys:
                c = cache.get(k)
                if c and c[0] == "":                       # cached miss
                    if REMATCH_MISSES or as_int(year_of[k]) >= threshold:
                        recheck.add(k)
        to_match = new | recheck
        reused = len(cur_keys) - len(to_match)
        print(f"  cache: {len(cur_keys)} current keys; reuse {reused}; "
              f"new {len(new)}; recheck-miss {len(recheck)} "
              f"(seqdb {'changed' if seqdb_changed else 'unchanged'})")

    # (re)match the needed keys, per subtype
    by_sub = {}
    for (sub, name, reass, passage) in to_match:
        by_sub.setdefault(sub, []).append((name, reass, passage))
    for sub, triples in by_sub.items():
        tag = SUBTYPE_TAG.get(sub)
        if tag is None:
            continue
        if USE_INDEX:                                     # single-process name index
            hits = _match_indexed(tag, sorted(triples))
        else:                                             # multiprocess filter_name scan
            chunks = [c for c in chunked(sorted(triples), NPROC * 4) if c]
            hits = {}
            with Pool(NPROC, initializer=_init, initargs=(tag,)) as pool:
                for d in pool.map(_work, chunks):
                    hits.update(d)
        for (name, reass, passage) in triples:            # record hits AND misses
            k = (sub, name, reass, passage)
            if (name, reass, passage) in hits:
                sid, aa, nucl = hits[(name, reass, passage)]
                cache[k] = (sid, str(len(aa)), str(nucl), aa)
            else:
                cache[k] = ("", "", "", "")
        print(f"  {sub}: matched {sum(1 for t in triples if t in hits)}/"
              f"{len(triples)} newly-checked keys")

    # persist cache (prune to current keys so it can't grow unbounded)
    with open(CACHE, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["subtype", "name", "reassortant", "passage",
                    "seq_id", "aa_length", "nuc_length", "aa"])
        for (sub, name, reass, passage) in sorted(cur_keys):
            v = cache.get((sub, name, reass, passage), ("", "", "", ""))
            w.writerow([sub, name, reass, passage, *v])
    with open(META, "w") as f:
        f.write(f"seqdb {cur_hash}\nmode {mode}\n")

    # derive match.csv (hits only) for every current antigen
    matched = 0
    with open(os.path.join(CSV_DIR, "match.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["ag_id", "subtype", "seq_id", "aa_length", "nuc_length", "aa"])
        for ag_id, key, _ in rows:
            v = cache.get(key)
            if v and v[0]:
                sub = key[0]
                w.writerow([ag_id, sub, v[0], v[1], v[2], v[3]])
                matched += 1
    print(f"TOTAL  {matched}/{len(rows)} antigens with sequence "
          f"({100*matched/max(len(rows),1):.1f}%)  ->  out/csv/match.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
