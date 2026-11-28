"""The U14e case table (scripts/u14e/vv/cases.py): every case is a
declared table tag and every ladder tag has a case (a row can always
be regenerated), each route's column covers exactly the cases that
run it, the selection follows the table's order and refuses an
unknown tag, the declared expectations are the committed records'
(the terms' eigenvalues, the ladder's addresses), and the minimal
family's tag names its channels and declares the ground term with its
configuration."""

import pytest

from vv import cases
from vv.table import Route, table_keys

pytestmark = pytest.mark.fast


def test_every_case_is_a_declared_tag_and_every_ladder_tag_has_a_case():
    tags = {k.split("|")[2] for k in table_keys()} - {"-"}
    assert set(cases.CASES) == tags and len(cases.CASES) == 24
    for tag, case in cases.CASES.items():
        assert case.tag == tag
        for method in ("casci", "casscf", "sa-casscf"):
            # the free routes run on every window
            assert f"{method}|{Route.FREE.value}|{tag}" in table_keys()
        # the pinned and the constrained columns carry only the cases
        # that run them
        pinned = f"sa-casscf|{Route.PINNED.value}|{tag}" in table_keys()
        radial = f"sa-casscf|{Route.RADIAL.value}|{tag}" in table_keys()
        assert pinned is case.pinned and radial is case.radial


def test_selection_follows_the_table_order():
    order = [c.tag for c in cases.select_cases()]
    assert order[:4] == ["(4,12)", "(4,15)", "(4,16)", "(6,13)"]
    assert order[4:9] == ["(4,12)@r1", "(4,15)@r1", "(4,16)@r1", "(4,21)@r1",
                          "(6,13)@r1"]
    assert order[9:] == list(cases.LADDER_MIN)       # the minimal family, by size
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


def test_a_minimal_tag_names_its_channels():
    # the suffix lists what stands beside the base [5f, 6d]: one
    # channel per token, two s virtuals for 2s
    assert cases.channels_of("(4,12)m") == ()
    assert cases.channels_of("(4,13)m+s") == ("8s",)
    assert cases.channels_of("(4,14)m+2s") == ("8s", "9s")
    assert cases.channels_of("(4,26)m+f+d+2s") == ("5f", "7d", "8s", "9s")
    for tag in cases.LADDER_MIN:
        width = 12 + sum(2 * cases.CHANNELS[c].l + 1
                         for c in cases.channels_of(tag))
        assert width == int(tag.split(",")[1].split(")")[0])


def test_a_malformed_minimal_tag_is_refused():
    with pytest.raises(ValueError, match="not a minimal-family tag"):
        cases.channels_of("(4,12)")
    with pytest.raises(ValueError, match=r"unknown channel token\(s\) \['q'\]"):
        cases.channels_of("(4,15)m+q")
    # the tag's size must be the channels' widths beside the 12-orbital base
    with pytest.raises(ValueError, match="make 15 active orbitals"):
        cases.channels_of("(4,16)m+p")


def test_the_minimal_family_declares_the_ground_term_and_its_configuration():
    # <L^2> and <S^2> cannot tell one configuration of the 5L term from
    # another (the frozen key window's is 5f^2 6d^2): the family gates
    # the electrons of the base's f and d blocks, and keeps its
    # composition only on the constrained route
    for tag in cases.LADDER_MIN:
        case = cases.CASES[tag]
        assert case.frozen == (cases.L2_TERM, cases.S2_QUINTET)
        assert case.ground_config == (3.0, 1.0)
        assert case.radial and case.radial_term == (cases.L2_TERM,
                                                    cases.S2_QUINTET)
        assert not case.pinned and not case.start_closed
        assert case.casscf_cycles == 50       # its single root runs to the cap
        assert case.sa_term is None            # the free SA converts s/p/d
    others = [c for t, c in cases.CASES.items() if t not in cases.LADDER_MIN]
    assert all(not c.radial and c.pinned and c.start_closed
               and c.ground_config is None and c.casscf_cycles is None
               for c in others)
