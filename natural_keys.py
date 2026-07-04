"""Content-based natural keys for antigens / sera / titer tables.

Positional hidb5 indices (`{sub}:{i}`) are NOT stable across hidb regenerations —
the array order shifts, so an id that means antigen X today can mean Y tomorrow.
That makes incremental append impossible and is why the titer build is full-rebuild.

These keys are derived from an entity's *identity content* instead, so the same
antigen/serum/table gets the same id across regenerations — the prerequisite for
incremental updates and the eventual Postgres system-of-record.

Key composition (validated against hidb5 for 0 titer-grain conflicts over 3.58M
cells; see the migration notes in README):

  antigen  {sub}:a:{h}  h = hash(name, reassortant, annotations, passage)
  serum    {sub}:s:{h}  h = hash(serum_id, name, reassortant, annotations, passage, species)
  table    {sub}:t:{h}  h = hash(lab, assay, rbc, date, virus, content-hash)

where the antigen/serum key is ae's canonical "designation" (name + reassortant +
annotations + passage), and the table content-hash is REORDER-INVARIANT: the sorted
set of (ag_key, sr_key, titer) triples, so a row/column reshuffle across a hidb
regeneration doesn't change the table id. The `{sub}:` prefix is retained so
subtype partitioning and `split_part(id, ':', 1)` keep working.

Shared by build_db.py (the producer) and build_clipped.py (which must reproduce the
same registered-entity keys) so the two can never drift.
"""
import hashlib
import json


def name_from(rec):
    """Reconstruct a strain name O/i/y, tolerating heterogeneous records.
    (Moved here from build_db so the key and the emitted `name` column agree.)"""
    iso = rec.get("i", "")
    if "/" in str(iso):           # already a full name (older serum records)
        return iso
    loc, yr = rec.get("O", ""), rec.get("y", "")
    return "/".join(p for p in (loc, iso, yr) if p)


def _annotations(rec):
    a = rec.get("a", [])
    if isinstance(a, list):
        return sorted(str(x) for x in a)
    return [str(a)] if a else []


def _h(*parts, n=16):
    return hashlib.md5(
        json.dumps(parts, sort_keys=True, ensure_ascii=False).encode("utf-8")
    ).hexdigest()[:n]


def antigen_key(sub, rec):
    return f"{sub}:a:{_h(name_from(rec), rec.get('R', ''), _annotations(rec), rec.get('P', ''))}"


def serum_key(sub, rec):
    return f"{sub}:s:{_h(rec.get('I', ''), name_from(rec), rec.get('R', ''), _annotations(rec), rec.get('P', ''), rec.get('s', ''))}"


def table_content_hash(rec, ag_keys, sr_keys):
    """Reorder-invariant hash of a table's registered (ag_key, sr_key, titer) cells."""
    ai, si, matrix = rec["a"], rec["s"], rec["t"]
    cells = []
    for r, row in enumerate(matrix):
        if r >= len(ai):
            continue
        ak = ag_keys[ai[r]]
        for c, cell in enumerate(row):
            if c >= len(si):
                continue
            s = str(cell).strip()
            if s and s != "*":
                cells.append((ak, sr_keys[si[c]], s))
    return _h(sorted(cells), n=12)


def table_key(sub, rec, content_hash):
    return f"{sub}:t:{_h(rec.get('l', ''), rec.get('A', ''), rec.get('r', ''), str(rec.get('D', '')), rec.get('V', ''), content_hash)}"
