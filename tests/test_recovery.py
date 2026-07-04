"""Invariants for the clipped-cell recovery (build_clipped.py -> load_duckdb.py).

These assert the provenance model holds in the built DB: every row is tagged
`hidb` or `ace_recovered`; recovered rows use the content-derived id namespace and
can't collide with hidb's; recovered titers stay FK-clean, attach to real tables,
and follow the same censoring convention. All skip cleanly if the DB isn't built
or if recovery wasn't run (ae_backend / whocc-tables absent) — so they never fail
a titers-only build.
"""
import math

import pytest


@pytest.fixture
def has_recovery(db):
    n = db.execute("SELECT count(*) FROM titer WHERE source='ace_recovered'").fetchone()[0]
    if n == 0:
        pytest.skip("no recovered data in this build (build_clipped not run)")
    return n


def test_source_domain(db):
    for tbl in ("antigen", "serum", "titer"):
        bad = db.execute(
            f"SELECT count(*) FROM {tbl} WHERE source NOT IN ('hidb','ace_recovered')"
        ).fetchone()[0]
        assert bad == 0, f"{tbl} has rows with an unexpected source"


def test_recovered_ids_namespaced(db, has_recovery):
    # recovered antigen/serum ids look like {sub}:ra:hash / {sub}:rs:hash and never
    # collide with hidb's positional {sub}:{int}
    bad_ag = db.execute(
        "SELECT count(*) FROM antigen WHERE source='ace_recovered' "
        "AND ag_id NOT SIMILAR TO '(h1|h3|b):ra:[0-9a-f]+'").fetchone()[0]
    assert bad_ag == 0
    bad_sr = db.execute(
        "SELECT count(*) FROM serum WHERE source='ace_recovered' "
        "AND sr_id NOT SIMILAR TO '(h1|h3|b):rs:[0-9a-f]+'").fetchone()[0]
    assert bad_sr == 0
    # no id appears under both provenances
    clash = db.execute(
        "SELECT count(*) FROM (SELECT ag_id FROM antigen GROUP BY ag_id "
        "HAVING count(DISTINCT source) > 1)").fetchone()[0]
    assert clash == 0


def test_recovered_titers_fk_clean(db, has_recovery):
    orphans = db.execute("""
        SELECT count(*) FROM titer t
        LEFT JOIN antigen a USING(ag_id)
        LEFT JOIN serum s USING(sr_id)
        LEFT JOIN titer_table tb USING(tab_id)
        WHERE t.source='ace_recovered'
          AND (a.ag_id IS NULL OR s.sr_id IS NULL OR tb.tab_id IS NULL)
    """).fetchone()[0]
    assert orphans == 0


def test_recovered_titers_attach_to_existing_tables(db, has_recovery):
    # recovery adds no new titer_table rows — recovered titers reuse hidb tab_ids
    new_tabs = db.execute(
        "SELECT count(DISTINCT tab_id) FROM titer WHERE source='ace_recovered' "
        "AND tab_id NOT IN (SELECT tab_id FROM titer_table)").fetchone()[0]
    assert new_tabs == 0


def test_recovered_titers_follow_thresholded_convention(db, has_recovery):
    bad = db.execute("""
        SELECT count(*) FROM titer WHERE source='ace_recovered' AND titer_kind<>'other'
        AND log_titer_thresholded <> CASE titer_kind
            WHEN 'lt' THEN log_titer-1 WHEN 'gt' THEN log_titer+1 ELSE log_titer END
    """).fetchone()[0]
    assert bad == 0


def test_recovered_volume_is_sane(db, has_recovery):
    # recovery is ~0.9% of titers at most (the clipped fraction); a much larger
    # number would signal double-counting of already-loaded cells
    rec = has_recovery
    total = db.execute("SELECT count(*) FROM titer").fetchone()[0]
    assert 0 < rec < 0.03 * total
    # every recovered antigen carries a non-empty name (we never emit a blank id row)
    blank = db.execute(
        "SELECT count(*) FROM antigen WHERE source='ace_recovered' AND (name IS NULL OR name='')"
    ).fetchone()[0]
    assert blank == 0


def test_titer_flat_exposes_provenance(db, has_recovery):
    cols = {c[0] for c in db.execute("DESCRIBE titer_flat").fetchall()}
    assert {"source", "antigen_source"} <= cols
    assert db.execute(
        "SELECT count(*) FROM titer_flat WHERE source='ace_recovered'").fetchone()[0] > 0
