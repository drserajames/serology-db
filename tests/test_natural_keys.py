"""Unit tests for natural_keys — the content-based id scheme (issue #5).

These pin the properties the whole migration relies on: keys are deterministic
(same identity -> same id, the basis for cross-regeneration stability), they
discriminate the identity-bearing fields (name, reassortant, annotations,
passage), and the table content hash is reorder-invariant (a row/column reshuffle
across a hidb regeneration must NOT change a table's id).
"""
import re

import natural_keys as nk


# minimal hidb-shaped records
def AG(**kw):
    base = {"O": "TESTLOC", "i": "1", "y": "1777", "P": "MDCK1", "R": "", "a": []}
    base.update(kw)
    return base


class TestNameFrom:
    def test_reconstructs_from_parts(self):
        assert nk.name_from(AG()) == "TESTLOC/1/1777"

    def test_full_name_when_slashed(self):
        assert nk.name_from({"i": "PLACE/9/1888"}) == "PLACE/9/1888"


class TestAntigenKey:
    def test_format(self):
        assert re.fullmatch(r"h3:a:[0-9a-f]{16}", nk.antigen_key("h3", AG()))

    def test_deterministic(self):
        assert nk.antigen_key("h3", AG()) == nk.antigen_key("h3", AG())

    def test_subtype_scoped(self):
        assert nk.antigen_key("h1", AG()) != nk.antigen_key("h3", AG())

    def test_discriminates_passage(self):
        assert nk.antigen_key("h3", AG(P="E5")) != nk.antigen_key("h3", AG(P="MDCK1"))

    def test_discriminates_reassortant(self):
        assert nk.antigen_key("h3", AG(R="X-187")) != nk.antigen_key("h3", AG(R=""))

    def test_discriminates_annotations(self):
        # the fix for the recovered-titer grain collision: annotations matter
        assert nk.antigen_key("h3", AG(a=["EGG"])) != nk.antigen_key("h3", AG(a=[]))

    def test_annotation_order_invariant(self):
        assert nk.antigen_key("h3", AG(a=["X", "Y"])) == nk.antigen_key("h3", AG(a=["Y", "X"]))


class TestSerumKey:
    def SR(self, **kw):
        base = {"I": "F001", "O": "TESTLOC", "i": "1", "y": "1777", "P": "", "R": "",
                "a": [], "s": "ferret"}
        base.update(kw)
        return base

    def test_format_and_determinism(self):
        k = nk.serum_key("b", self.SR())
        assert re.fullmatch(r"b:s:[0-9a-f]{16}", k)
        assert k == nk.serum_key("b", self.SR())

    def test_discriminates_serum_id_and_species(self):
        assert nk.serum_key("b", self.SR(I="F002")) != nk.serum_key("b", self.SR(I="F001"))
        assert nk.serum_key("b", self.SR(s="cell")) != nk.serum_key("b", self.SR(s="ferret"))


class TestTableKey:
    def _tab(self, matrix, a, s, **meta):
        base = {"l": "LAB", "A": "HI", "r": "turkey", "D": "18880301", "V": "SEASON",
                "a": a, "s": s, "t": matrix}
        base.update(meta)
        return base

    def test_reorder_invariant_content_hash(self):
        ag_keys = ["ka", "kb", "kc"]
        sr_keys = ["s0", "s1"]
        # same biological content, rows + columns permuted
        t1 = self._tab([["40", "80"], ["160", "320"]], [0, 1], [0, 1])
        t2 = self._tab([["320", "160"], ["80", "40"]], [1, 0], [1, 0])
        h1 = nk.table_content_hash(t1, ag_keys, sr_keys)
        h2 = nk.table_content_hash(t2, ag_keys, sr_keys)
        assert h1 == h2

    def test_content_hash_changes_with_titer(self):
        ag_keys, sr_keys = ["ka"], ["s0"]
        a = self._tab([["40"]], [0], [0])
        b = self._tab([["80"]], [0], [0])
        assert nk.table_content_hash(a, ag_keys, sr_keys) != nk.table_content_hash(b, ag_keys, sr_keys)

    def test_table_key_format_and_metadata_sensitivity(self):
        t = self._tab([["40"]], [0], [0])
        k = nk.table_key("h1", t, "abc123")
        assert re.fullmatch(r"h1:t:[0-9a-f]{16}", k)
        assert nk.table_key("h1", dict(t, l="OTHER"), "abc123") != k
