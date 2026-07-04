"""Data-quality invariants asserted against the built serology.duckdb.

These are the "is the DB sane?" checks — the properties that must hold after any
refresh for the store to be trustworthy. They run read-only against the live DB
and skip cleanly if it hasn't been built (see the `db` fixture). Thresholds are
deliberately loose bands, not exact counts, so they survive a normal data refresh
but catch a structural regression (a broken join, a dropped subtype, a parsing
change that silently loses rows).
"""
import math


# --- referential integrity: no orphan facts --------------------------------
def test_no_orphan_titer_antigen(db):
    assert db.execute(
        "SELECT count(*) FROM titer t LEFT JOIN antigen a USING(ag_id) "
        "WHERE a.ag_id IS NULL").fetchone()[0] == 0


def test_no_orphan_titer_serum(db):
    assert db.execute(
        "SELECT count(*) FROM titer t LEFT JOIN serum s USING(sr_id) "
        "WHERE s.sr_id IS NULL").fetchone()[0] == 0


def test_no_orphan_titer_table(db):
    assert db.execute(
        "SELECT count(*) FROM titer t LEFT JOIN titer_table tb USING(tab_id) "
        "WHERE tb.tab_id IS NULL").fetchone()[0] == 0


def test_no_orphan_sequence_antigen(db):
    assert db.execute(
        "SELECT count(*) FROM sequence sq LEFT JOIN antigen a USING(ag_id) "
        "WHERE a.ag_id IS NULL").fetchone()[0] == 0


# --- content-based natural keys (issue #5) ---------------------------------
def test_ids_are_natural_keys(db):
    # hidb rows carry {sub}:a|s:{hash} keys and a positional hidb_id; recovered
    # rows carry {sub}:ra|rs:{hash} and NULL hidb_id.
    bad_ag = db.execute(
        "SELECT count(*) FROM antigen WHERE source='hidb' "
        "AND ag_id NOT SIMILAR TO '(h1|h3|b):a:[0-9a-f]+'").fetchone()[0]
    assert bad_ag == 0
    bad_tab = db.execute(
        "SELECT count(*) FROM titer_table "
        "WHERE tab_id NOT SIMILAR TO '(h1|h3|b):t:[0-9a-f]+'").fetchone()[0]
    assert bad_tab == 0
    assert db.execute("SELECT count(*) FROM antigen WHERE source='hidb' "
                      "AND hidb_id IS NULL").fetchone()[0] == 0


def test_titer_grain_unique(db):
    # the whole point of natural keys: (table, antigen, serum) is a real grain.
    dup = db.execute(
        "SELECT count(*) FROM (SELECT tab_id, ag_id, sr_id FROM titer "
        "GROUP BY ALL HAVING count(*) > 1)").fetchone()[0]
    assert dup == 0


# --- primary-key uniqueness ------------------------------------------------
def test_unique_primary_keys(db):
    for table, key in (("antigen", "ag_id"), ("serum", "sr_id"),
                       ("titer_table", "tab_id"), ("sequence", "ag_id"),
                       ("location", "location")):
        dupes = db.execute(
            f"SELECT count(*) FROM (SELECT {key} FROM {table} "
            f"GROUP BY {key} HAVING count(*) > 1)").fetchone()[0]
        assert dupes == 0, f"{table}.{key} has {dupes} duplicated keys"


# --- volume bands: catch a dropped subtype / halved load -------------------
def test_row_count_bands(db):
    counts = {t: db.execute(f"SELECT count(*) FROM {t}").fetchone()[0]
              for t in ("antigen", "serum", "titer_table", "titer",
                        "sequence", "location")}
    # loose lower/upper bounds around the known ~2026 baseline
    assert 200_000 <= counts["antigen"] <= 500_000
    assert 3_000 <= counts["serum"] <= 20_000
    assert 5_000 <= counts["titer_table"] <= 30_000
    assert 2_500_000 <= counts["titer"] <= 8_000_000
    assert 50_000 <= counts["sequence"] <= counts["antigen"]
    assert 2_000 <= counts["location"] <= 20_000


def test_all_three_subtypes_present(db):
    subs = {r[0] for r in db.execute(
        "SELECT DISTINCT subtype FROM antigen").fetchall()}
    assert {"h1", "h3", "b"} <= subs


# --- titer partition + value/log sanity ------------------------------------
def test_titer_kind_partition(db):
    kinds = dict(db.execute(
        "SELECT titer_kind, count(*) FROM titer GROUP BY ALL").fetchall())
    assert set(kinds) <= {"num", "lt", "gt", "other"}
    total = db.execute("SELECT count(*) FROM titer").fetchone()[0]
    assert sum(kinds.values()) == total
    # exact readings must dominate
    assert kinds.get("num", 0) > 0.5 * total


def test_log_titer_consistency(db):
    # every non-'other' titer has a value and a log; 'other' has neither
    bad = db.execute(
        "SELECT count(*) FROM titer WHERE titer_kind <> 'other' "
        "AND (titer_value IS NULL OR log_titer IS NULL)").fetchone()[0]
    assert bad == 0
    other_with_val = db.execute(
        "SELECT count(*) FROM titer WHERE titer_kind = 'other' "
        "AND titer_value IS NOT NULL").fetchone()[0]
    assert other_with_val == 0


def test_log_titer_matches_value(db):
    # log_titer == log2(value/10) must hold for a large sample
    rows = db.execute(
        "SELECT titer_value, log_titer FROM titer "
        "WHERE titer_value IS NOT NULL USING SAMPLE 5000 ROWS").fetchall()
    assert rows
    for val, lg in rows:
        assert lg == math.log2(val / 10.0)


def test_log_titer_thresholded_convention(db):
    # ae Titer::logged_with_thresholded(): <N -> log-1 (N/2), >N -> log+1 (N*2),
    # regular unchanged. No non-'other' row may violate this.
    bad = db.execute("""
        SELECT count(*) FROM titer WHERE titer_kind <> 'other' AND
          log_titer_thresholded <> CASE titer_kind
            WHEN 'lt' THEN log_titer - 1
            WHEN 'gt' THEN log_titer + 1
            ELSE log_titer END
    """).fetchone()[0]
    assert bad == 0
    # and it must differ from face value exactly on the censored rows
    diff = db.execute(
        "SELECT count(*) FROM titer WHERE log_titer_thresholded <> log_titer"
    ).fetchone()[0]
    censored = db.execute(
        "SELECT count(*) FROM titer WHERE titer_kind IN ('lt','gt')"
    ).fetchone()[0]
    assert diff == censored


def test_log_titer_in_plausible_range(db):
    lo, hi = db.execute(
        "SELECT min(log_titer), max(log_titer) FROM titer").fetchone()
    assert -4.0 <= lo <= 0.0
    assert 6.0 <= hi <= 20.0


# --- resolution / coverage rates -------------------------------------------
def test_location_resolution_rate(db):
    resolved, total = db.execute(
        "SELECT count(*) FILTER(WHERE country <> ''), count(*) "
        "FROM location").fetchone()
    assert resolved / total >= 0.98


def test_antigen_location_row_resolution(db):
    # share of antigen rows whose raw location resolves to a country
    resolved, total = db.execute(
        "SELECT count(*) FILTER(WHERE loc.country IS NOT NULL "
        "AND loc.country <> ''), count(*) "
        "FROM antigen a LEFT JOIN location loc ON a.location = loc.location "
        "WHERE a.location <> ''").fetchone()
    assert resolved / total >= 0.99


def test_sequence_coverage_band(db):
    ag = db.execute("SELECT count(*) FROM antigen").fetchone()[0]
    seq = db.execute("SELECT count(*) FROM sequence").fetchone()[0]
    # documented ~34-44%; band leaves headroom either side
    assert 0.25 <= seq / ag <= 0.60


# --- views resolve and join correctly --------------------------------------
def test_titer_flat_view_matches_titer_count(db):
    # titer_flat inner-joins titer to antigen/serum/table; with FK integrity it
    # must preserve the titer row count exactly (LEFT joins to seq/loc don't fan
    # out because ag_id/location are unique in those tables).
    tf = db.execute("SELECT count(*) FROM titer_flat").fetchone()[0]
    ti = db.execute("SELECT count(*) FROM titer").fetchone()[0]
    assert tf == ti


def test_antigen_sequence_view_matches_antigen_count(db):
    av = db.execute("SELECT count(*) FROM antigen_sequence").fetchone()[0]
    a = db.execute("SELECT count(*) FROM antigen").fetchone()[0]
    assert av == a
