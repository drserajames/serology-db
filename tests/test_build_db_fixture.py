"""End-to-end regression test of build_db over a synthetic hidb5.

Runs the real build_db.main() against a tiny, fully-fake hidb5.h1.json.xz (see
conftest.synthetic_hidb5) and asserts the emitted CSVs exactly. This locks down
the reshape logic that Phase C (clipped-cell reconciliation) will touch: if a
change alters how many titers load or how clipped cells are counted, this fails.
"""
import csv
import os

from conftest import EXPECTED


def _read_csv(csv_dir, name):
    with open(os.path.join(csv_dir, name)) as f:
        return list(csv.DictReader(f))


def test_row_counts_and_clipping(fixture_env, capsys):
    build_db = fixture_env["module"]
    csv_dir = str(fixture_env["csv_dir"])

    n_ti = build_db.main()
    out = capsys.readouterr().out

    antigens = _read_csv(csv_dir, "antigen.csv")
    sera = _read_csv(csv_dir, "serum.csv")
    tables = _read_csv(csv_dir, "titer_table.csv")
    titers = _read_csv(csv_dir, "titer.csv")

    assert len(antigens) == EXPECTED["n_antigens"]
    assert len(sera) == EXPECTED["n_sera"]
    assert len(tables) == EXPECTED["n_tables"]
    assert len(titers) == EXPECTED["n_titers"]
    assert n_ti == EXPECTED["n_titers"]
    # clipped-cell tally is reported to stdout, not returned
    assert f"clipped_cells={EXPECTED['n_clipped']}" in out


def test_ids_are_subtype_prefixed(fixture_env):
    build_db = fixture_env["module"]
    csv_dir = str(fixture_env["csv_dir"])
    build_db.main()

    antigens = _read_csv(csv_dir, "antigen.csv")
    assert {a["ag_id"] for a in antigens} == {"h1:0", "h1:1", "h1:2"}
    titers = _read_csv(csv_dir, "titer.csv")
    for t in titers:
        assert t["tab_id"] == "h1:0"
        assert t["ag_id"].startswith("h1:")
        assert t["sr_id"].startswith("h1:")


def test_titer_kinds_and_values(fixture_env):
    build_db = fixture_env["module"]
    csv_dir = str(fixture_env["csv_dir"])
    build_db.main()

    titers = _read_csv(csv_dir, "titer.csv")
    by_raw = {t["titer_raw"]: t for t in titers}
    # exactly the three non-clipped, non-missing cells survive
    assert set(by_raw) == {"1280", "<10", "640"}
    assert by_raw["1280"]["titer_kind"] == "num"
    assert by_raw["<10"]["titer_kind"] == "lt"
    assert by_raw["640"]["titer_kind"] == "num"
    # the '*' cell at (1,0) must be absent
    assert all(t["titer_raw"] != "*" for t in titers)


def test_table_date_parsing(fixture_env):
    build_db = fixture_env["module"]
    csv_dir = str(fixture_env["csv_dir"])
    build_db.main()

    tables = _read_csv(csv_dir, "titer_table.csv")
    assert tables[0]["table_date"] == "1888-03-01"
    assert tables[0]["assay"] == "HI"
    assert tables[0]["lab"] == "TESTLAB"


def test_serum_name_reconstruction(fixture_env):
    build_db = fixture_env["module"]
    csv_dir = str(fixture_env["csv_dir"])
    build_db.main()

    sera = _read_csv(csv_dir, "serum.csv")
    names = {s["serum_id"]: s["name"] for s in sera}
    assert names["F0001"] == "TESTLOC/1/1777"        # reconstructed from O/i/y
    assert names["F0002"] == "SAMPLETOWN/9/1888"     # taken verbatim from slashed i
