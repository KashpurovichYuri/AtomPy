"""The U14e case table (scripts/u14e/vv/cases.py): every case is a
declared table tag and every ladder tag has a case (a row can always
be regenerated), the selection follows the table's order and refuses
an unknown tag, and the declared expectations are the committed
records' (the terms' eigenvalues, the ladder's addresses)."""

import pytest

from vv import cases
from vv.table import Route, table_keys

pytestmark = pytest.mark.fast


def test_every_case_is_a_declared_tag_and_every_ladder_tag_has_a_case():
    tags = {k.split("|")[2] for k in table_keys()} - {"-"}
    assert set(cases.CASES) == tags
    for tag, case in cases.CASES.items():
        assert case.tag == tag
        for method in ("casci", "casscf", "sa-casscf"):
            assert f"{method}|{Route.FREE.value}|{tag}" in table_keys()
        assert f"sa-casscf|{Route.PINNED.value}|{tag}" in table_keys()


def test_selection_follows_the_table_order():
    order = [c.tag for c in cases.select_cases()]
    assert order[:4] == ["(4,12)", "(4,15)", "(4,16)", "(6,13)"]
    assert order[4:] == ["(4,12)@r1", "(4,15)@r1", "(4,16)@r1", "(4,21)@r1",
                         "(6,13)@r1"]
    assert [c.tag for c in cases.select_cases(["(6,13)", "(4,12)"])] == \
        ["(6,13)", "(4,12)"]
    assert [c.tag for c in cases.select_cases(["(4,12)", "(4,12)"])] == ["(4,12)"]
    with pytest.raises(ValueError, match="unknown case tag"):
        cases.select_cases(["(4,18)"])


def test_the_declared_terms_and_gates():
    assert cases.TERM == (8, "u") and cases.NROOTS == 17
    assert (cases.L2_TERM, cases.L2_L9) == (72.0, 90.0)
    assert (cases.S2_QUINTET, cases.S2_SEPTET) == (6.0, 12.0)
    gate = cases.CASES["(4,12)@r1"]
    assert gate.pinned_frozen == gate.pinned_sa == cases.L2_TERM
    assert gate.frozen == (cases.L2_TERM, cases.S2_QUINTET)
    six = cases.CASES["(6,13)"]
    assert six.frozen == (12.0, cases.S2_SEPTET) and six.pinned_class == (16, 1)
    assert cases.CASES["(6,13)@r1"].s_pair
    assert cases.CASES["(6,13)@r1"].frozen_below == 400.0
    assert cases.CASES["(4,12)"].anchor == "u14e_10f13d_4-0_sa-casscf_ao_-_pur_avg"
    # the pin-on-pur gate is declared for exactly one case
    assert [t for t, c in cases.CASES.items() if c.pinned_sa is not None] == ["(4,12)@r1"]


def test_the_addresses_are_unique():
    assert len(set(cases.ADDR.values())) == len(cases.ADDR)
