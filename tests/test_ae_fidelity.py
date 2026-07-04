"""Scientific-fidelity check: our titer parsing vs ae's canonical Titer class.

This reproduces the reference computation on a *sample of real charts* and asserts
that, for every distinct raw titer string encountered:

  * build_db.parse_titer's log        == ae Titer.logged()               (face value)
  * the thresholded convention we apply == ae Titer.logged_with_thresholded()
        (i.e. plain log for regulars, log-1 for '<N', log+1 for '>N')
  * censoring classification agrees   ('<' -> lt, '>' -> gt, '*' -> omitted)

It is deliberately NOT hermetic: it needs the built ae_backend extension and the
read-only whocc-tables charts. Both are skipped-with-reason when absent, so the
core suite stays dependency-free and portable. When they ARE present this is the
authoritative guard that our tidy log-titers match the C++ toolkit byte-for-byte.
"""
import glob
import math
import os

import pytest

AE_BUILD = os.environ.get(
    "AE_BUILD",
    os.path.expanduser("~/AC/eu/ae/build"),
)
WHOCC = os.environ.get(
    "WHOCC_TABLES",
    os.path.expanduser("~/AC/eu/whocc-tables"),
)


@pytest.fixture(scope="module")
def ae_titer():
    import sys
    sys.path.insert(0, AE_BUILD)
    try:
        import ae_backend  # noqa: F401
    except Exception as err:  # pragma: no cover - environment dependent
        pytest.skip(f"ae_backend unavailable: {err}")
    return ae_backend


@pytest.fixture(scope="module")
def sample_charts():
    charts = sorted(glob.glob(os.path.join(WHOCC, "**", "*.ace"), recursive=True))
    if not charts:
        pytest.skip(f"no .ace charts under {WHOCC}")
    # a small, deterministic spread across subtype/assay dirs — every 500th chart
    return charts[::500][:12] or charts[:4]


def _thresholded(parsed):
    """Apply the load_duckdb convention to a parse_titer result -> logged_with_thresholded."""
    kind, _val, lg = parsed
    if lg is None:
        return None
    if kind == "lt":
        return lg - 1
    if kind == "gt":
        return lg + 1
    return lg


def test_parse_titer_matches_ae_logged(ae_titer, sample_charts):
    import build_db

    c3 = ae_titer.chart_v3
    checked = 0
    censored_seen = 0
    for path in sample_charts:
        chart = c3.Chart(path)
        titers = chart.titers()
        na, ns = chart.number_of_antigens(), chart.number_of_sera()
        for ai in range(na):
            for si in range(ns):
                t = titers.titer(ai, si)
                raw = str(t)
                if t.is_dont_care():
                    # '*' / missing: ae throws for logged(); we return None (omit)
                    assert build_db.parse_titer(raw) is None
                    continue
                parsed = build_db.parse_titer(raw)
                assert parsed is not None, f"we dropped a non-missing titer {raw!r}"
                kind, _val, lg = parsed
                # face-value log matches ae logged()
                assert lg == pytest.approx(t.logged()), \
                    f"{raw!r}: our log {lg} != ae logged {t.logged()}"
                # thresholded convention matches ae logged_with_thresholded()
                assert _thresholded(parsed) == pytest.approx(t.logged_with_thresholded()), \
                    f"{raw!r}: our thresholded {_thresholded(parsed)} != ae {t.logged_with_thresholded()}"
                # censoring classification agrees
                assert (kind == "lt") == t.is_less_than(), f"{raw!r} lt mismatch"
                assert (kind == "gt") == t.is_more_than(), f"{raw!r} gt mismatch"
                if kind in ("lt", "gt"):
                    censored_seen += 1
                checked += 1
    assert checked > 0, "no titers checked — sample charts were empty"
    # the whole point is the censored path; make sure we actually exercised it
    assert censored_seen > 0, "no censored titers in the sample — widen the sample"
