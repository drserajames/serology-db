"""Shared pytest fixtures for the serology-db test suite.

Two kinds of tests live here:
  * unit / fixture tests — hermetic, build a tiny synthetic hidb5 through the real
    pipeline code in a tmp dir. No real data, no built DB required. Always run.
  * invariant tests — assert data-quality properties of the *built* serology.duckdb
    (row counts, FK integrity, censoring, resolution rates). Skipped with a clear
    message if the DB hasn't been built yet.

The synthetic fixtures contain only obviously-fake values and are written to a
tmp dir, never to the repo — no WHO data is ever committed.
"""
import json
import lzma
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
# The build scripts live at the repo root; make them importable.
if REPO not in sys.path:
    sys.path.insert(0, REPO)
ACMACS_DATA = os.environ.get(
    "ACMACS_DATA", os.path.normpath(os.path.join(REPO, os.pardir, "acmacs-data"))
)
SERO_OUT = os.environ.get("SERO_OUT", os.path.join(ACMACS_DATA, "serology-db"))
DB_PATH = os.path.join(SERO_OUT, "serology.duckdb")


# ---------------------------------------------------------------------------
# Synthetic hidb5 fixture — a minimal but structurally faithful subtype file.
# ---------------------------------------------------------------------------
def synthetic_hidb5():
    """A tiny, fully-fake hidb5 subtype dict exercising the pipeline's edge cases.

    Designed so the expected build_db output is known exactly:

    Antigens (3), sera (2), one table indexing antigens [0,1] and sera [0,1] but
    carrying a 3x3 matrix — so the trailing row/column are 'clipped' (titer cells
    with no registered antigen/serum index). The matrix is chosen so:

        loaded titers      = 3   (num, lt, and one more num)
        clipped real cells = 4   (2 in the extra column, 2 in the extra row)
        '*' cells          = omitted everywhere (never loaded, never clipped)
    """
    # Locations are invented and years are deliberately outside the 19xx/20xx
    # range so these fake strain names cannot match the pre-commit guardrail's
    # real-data heuristic (LOCATION/iso/19xx|20xx). The shape is still faithful.
    return {
        "a": [
            {"O": "TESTLOC", "i": "1", "y": "1777", "V": "A(H1N1)",
             "L": "", "P": "MDCK1", "R": "", "D": ["1777-01-15"]},
            {"O": "WI", "i": "2", "y": "1888", "V": "A(H1N1)",
             "L": "", "P": "E3", "R": "", "D": ["1888-02-20"]},
            {"O": "OTHERTOWN", "i": "3", "y": "1650", "V": "A(H1N1)",
             "L": "", "P": "SIAT1", "R": "", "D": []},
        ],
        "s": [
            # serum whose name is reconstructed from O/i/y (no slash in i)
            {"I": "F0001", "O": "TESTLOC", "i": "1", "y": "1777",
             "V": "A(H1N1)", "L": "", "P": "", "s": "ferret"},
            # serum whose `i` already holds a full slashed name -> used verbatim
            {"I": "F0002", "O": "", "i": "SAMPLETOWN/9/1888", "y": "",
             "V": "A(H1N1)", "L": "", "P": "", "s": "ferret"},
        ],
        "t": [
            {"A": "HI", "r": "turkey", "l": "TESTLAB", "V": "SEASON",
             "D": "18880301",
             "a": [0, 1], "s": [0, 1],
             "t": [
                 ["1280", "<10", "40"],      # row0: (num, lt) loaded; "40" clipped col
                 ["*", "640", ">2560"],       # row1: "*" omit, 640 loaded, ">2560" clipped col
                 ["80", "*", "160"],          # row2: whole row clipped -> "80","160" count, "*" not
             ]},
        ],
    }


EXPECTED = {
    "n_antigens": 3,
    "n_sera": 2,
    "n_tables": 1,
    "n_titers": 3,
    "n_clipped": 4,
}


@pytest.fixture
def fixture_env(tmp_path, monkeypatch):
    """Run build_db against a synthetic hidb5.h1.json.xz in an isolated tmp dir.

    Yields the tmp acmacs-data dir and csv output dir. Imports build_db and points
    its module-level path globals at the tmp locations, then writes a real xz file
    so load_json_xz (incl. the \\U CJK-escape normalisation) is exercised too.
    """
    import build_db

    data_dir = tmp_path / "acmacs-data"
    csv_dir = tmp_path / "out" / "csv"
    data_dir.mkdir(parents=True)
    csv_dir.mkdir(parents=True)

    payload = json.dumps(synthetic_hidb5())
    with lzma.open(data_dir / "hidb5.h1.json.xz", "wb") as f:
        f.write(payload.encode("utf-8"))

    monkeypatch.setattr(build_db, "ACMACS_DATA", str(data_dir))
    monkeypatch.setattr(build_db, "CSV_DIR", str(csv_dir))
    monkeypatch.setattr(build_db, "SUBTYPES", ["h1"])
    return {"module": build_db, "data_dir": data_dir, "csv_dir": csv_dir}


# ---------------------------------------------------------------------------
# Live built-DB connection (read-only). Skips cleanly if the DB isn't built.
# ---------------------------------------------------------------------------
@pytest.fixture(scope="session")
def db():
    if not os.path.exists(DB_PATH):
        pytest.skip(f"built DB not found at {DB_PATH}; run ./refresh.sh first")
    duckdb = pytest.importorskip("duckdb")
    con = duckdb.connect(DB_PATH, read_only=True)
    yield con
    con.close()
