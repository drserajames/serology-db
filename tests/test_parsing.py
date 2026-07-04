"""Unit tests for build_db's pure transformation functions.

These are the core of the tidy-ification: titer parsing (kind/value/log),
strain-name reconstruction, and the list-flattening helper. They have no I/O and
no data dependency, so they run everywhere and pin the exact semantics the rest
of the pipeline (and the DB schema) relies on.
"""
import math

import pytest

import build_db


class TestParseTiter:
    def test_plain_numeric(self):
        assert build_db.parse_titer("1280") == ("num", 1280, math.log2(128))

    def test_tilde_is_numeric(self):
        # '~' is an approximate/estimated reading — treated as an exact num value
        kind, val, lg = build_db.parse_titer("~80")
        assert (kind, val) == ("num", 80)
        assert lg == pytest.approx(3.0)

    def test_left_censored(self):
        # NOTE: current convention stores the threshold value itself (value=10),
        # so log_titer==0. Phase B revisits whether <10 should map to log -1.
        assert build_db.parse_titer("<10") == ("lt", 10, 0.0)

    def test_right_censored(self):
        assert build_db.parse_titer(">2560") == ("gt", 2560, math.log2(256))

    @pytest.mark.parametrize("raw", ["*", "", "   ", None])
    def test_missing_is_none(self, raw):
        assert build_db.parse_titer(raw) is None

    @pytest.mark.parametrize("raw", ["abc", "ND", "QNS"])
    def test_unparseable_is_other(self, raw):
        assert build_db.parse_titer(raw) == ("other", None, None)

    def test_whitespace_is_stripped(self):
        assert build_db.parse_titer("  40 ") == ("num", 40, 2.0)

    def test_low_titer_gives_negative_log(self):
        kind, val, lg = build_db.parse_titer("5")
        assert (kind, val, lg) == ("num", 5, -1.0)

    def test_numeric_types(self):
        # hidb matrices sometimes carry ints, not strings
        assert build_db.parse_titer(160) == ("num", 160, 4.0)


class TestNameFrom:
    # Fake locations + non-19xx/20xx years keep these clear of the pre-commit
    # guardrail (which flags real LOCATION/iso/19xx|20xx strain names).
    def test_reconstructs_from_parts(self):
        rec = {"O": "TESTLOC", "i": "1", "y": "1777"}
        assert build_db.name_from(rec) == "TESTLOC/1/1777"

    def test_uses_full_name_when_i_is_slashed(self):
        rec = {"O": "IGNORED", "i": "SAMPLETOWN/9/1888", "y": "1650"}
        assert build_db.name_from(rec) == "SAMPLETOWN/9/1888"

    def test_skips_empty_parts(self):
        assert build_db.name_from({"O": "OTHERTOWN", "i": "", "y": "1650"}) == "OTHERTOWN/1650"

    def test_all_empty(self):
        assert build_db.name_from({"O": "", "i": "", "y": ""}) == ""


class TestFirst:
    def test_first_of_list(self):
        assert build_db.first(["2019-01-15", "x"]) == "2019-01-15"

    def test_empty_list_is_none(self):
        assert build_db.first([]) is None

    def test_scalar_passthrough(self):
        assert build_db.first("scalar") == "scalar"

    def test_falsey_scalar_is_none(self):
        assert build_db.first("") is None
