"""Fidelity test for build_sequences' opt-in --index fast path.

Asserts that `_match_indexed` returns exactly what the default `filter_name`
scan does, on a real sample of antigen keys. Like test_ae_fidelity this is not
hermetic — it needs the built ae_backend seqdb and a built antigen.csv — so it
skips cleanly when either is absent. When they're present it's the guard that the
name-index path can never silently diverge from the canonical matcher.
"""
import csv
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
AE_BUILD = os.environ.get("AE_BUILD", os.path.expanduser("~/AC/eu/ae/build"))
ACMACS = os.environ.get(
    "ACMACS_DATA", os.path.normpath(os.path.join(REPO, os.pardir, "acmacs-data")))
SERO_OUT = os.environ.get("SERO_OUT", os.path.join(ACMACS, "serology-db"))
ANTIGEN = os.path.join(os.environ.get("SERO_CSV_DIR", os.path.join(SERO_OUT, "csv")),
                       "antigen.csv")
csv.field_size_limit(1 << 24)


@pytest.fixture(scope="module")
def ae():
    sys.path.insert(0, AE_BUILD)
    try:
        import ae_backend  # noqa: F401
    except Exception as err:  # pragma: no cover - env dependent
        pytest.skip(f"ae_backend unavailable: {err}")
    return ae_backend


def _sample_keys(subtype, tag, limit=500):
    if not os.path.exists(ANTIGEN):
        pytest.skip(f"no antigen.csv at {ANTIGEN}; run the pipeline first")
    seen, keys = set(), []
    with open(ANTIGEN) as f:
        for r in csv.DictReader(f):
            if r["subtype"] != subtype:
                continue
            k = (r["name"], r["reassortant"], "")
            if k not in seen:
                seen.add(k)
                keys.append(k)
    # deterministic spread across the subtype's key space
    step = max(1, len(keys) // limit)
    return keys[::step][:limit]


# h3 only: the largest subtype and the one with the most multi-candidate names
# (the filter_name-fallback path), so it exercises every branch. Keeping it to one
# subtype avoids loading three seqdbs twice — the full three-subtype match.csv was
# separately verified byte-identical to the slow path during development.
def test_index_path_matches_filter_name(ae):
    import build_sequences as bs

    subtype, tag = "h3", "A(H3N2)"
    triples = _sample_keys(subtype, tag, limit=400)
    assert triples, "no keys sampled"

    got = bs._match_indexed(tag, triples)

    # reference: the default per-key select_all().filter_name path
    sdb = ae.seqdb.for_subtype(tag)
    ref = {}
    for (name, reass, passage) in triples:
        sel = sdb.select_all().filter_name(name=name, reassortant=reass, passage=passage)
        if len(sel):
            r = sel[0]
            ref[(name, reass, passage)] = (r.seq_id(), str(r.aa), len(str(r.nuc)))

    assert got == ref            # identical hits, seq_ids, aa, and nuc lengths
    assert ref, "sample produced no hits — widen it"
