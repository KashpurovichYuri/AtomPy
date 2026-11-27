"""The V&V table module (scripts/u14e/vv/table.py) on a temporary store
(the committed one is never touched): the store accepts only the
declared keys, the table renders one block per method class with a rule
row between blocks and the plain ladder before the @r1 rows, a
withdrawn row renders dashes in the withdrawn cells only, and the
bookkeeping E_tot = E_method + E(2) is checked at render."""

import json

import pytest

from vv import table

pytestmark = pytest.mark.fast


@pytest.fixture
def store(tmp_path, monkeypatch):
    """The row store and the table redirected to tmp_path; the ROHF
    energy (the zero of every dE) fixed."""
    monkeypatch.setattr(table, "ROWS", tmp_path / "rows.json")
    monkeypatch.setattr(table, "TABLE", tmp_path / "table.md")
    monkeypatch.setattr(table, "rohf_energy", lambda: -50.0)
    return tmp_path


def _row(method, route, cas, e_method=-50.4, e2_meh=-300.0, l2="72.0000"):
    # a reference row (no PT2 by definition) totals to its own energy
    e_tot = (e_method if route == "-"
             else None if e2_meh is None else e_method + e2_meh * 1e-3)
    return table.std_row(method, route, cas, e_method, e2_meh, e_tot, l2, "6.0000")


def _body(tmp_path):
    """The data rows and the rule rows (the title, the blank line, the
    header and its rule are skipped)."""
    return (tmp_path / "table.md").read_text().splitlines()[4:]


def _cells(line):
    return [c.strip() for c in line.split("|")[1:-1]]


def test_drop_row_refuses_an_undeclared_key(store):
    for row in (_row("casci", "free", "(4,12)"),          # the old spelling
                _row("casci", "free@pur", "(4,99)"),      # off the ladder
                _row("sa-casscf", "pinned@pur", "(4,18)"),
                _row("casci", "pinned@pur", "(4,12)")):   # not a kept route
        with pytest.raises(ValueError, match="outside the declared table"):
            table.drop_row(row)
    assert not (store / "rows.json").exists()             # nothing was written


def test_render_refuses_an_undeclared_key_in_the_store(store):
    key = "casci|free|(4,12)"
    (store / "rows.json").write_text(json.dumps({key: _row("casci", "free@pur", "(4,12)")}))
    with pytest.raises(ValueError, match=r"casci\|free\|\(4,12\)"):
        table.render_table()


def test_one_block_per_method_class_with_a_rule_between(store):
    for m in ("sa-casscf", "casci", "casscf"):             # dropped out of order
        table.drop_row(_row(m, "free@pur", "(4,12)"))
    table.drop_row(_row("sa-casscf", "pinned@pur", "(4,12)"))
    table.drop_row(_row("rohf", "-", "-", e2_meh=None))
    body = _body(store)
    methods = [_cells(r)[0] for r in body if "===" not in r]
    assert methods == ["rohf", "casci (free@pur)", "casscf (free@pur)",
                       "sa-casscf (free@pur)", "sa-casscf (pinned@pur)"]
    assert [i for i, r in enumerate(body) if "===" in r] == [1, 3, 5]


def test_within_a_block_the_plain_ladder_precedes_the_r1_rows(store):
    for cas in ("(6,13)@r1", "(4,16)", "(4,12)@r1", "(4,15)"):
        table.drop_row(_row("casci", "free@pur", cas))
    assert [_cells(r)[1] for r in _body(store)] == [
        "(4,15)", "(4,16)", "(4,12)@r1", "(6,13)@r1"]   # size-ordered


def test_a_withdrawn_row_renders_dashes_in_the_withdrawn_cells_only(store):
    row = _row("sa-casscf", "pinned@pur", "(6,13)")
    for field in ("l2", "e2_meh", "e2_ev", "e_tot"):
        row[field] = None
    table.drop_row(row)
    c = _cells(_body(store)[0])
    assert (c[2], c[7], c[8], c[9]) == ("—", "—", "—", "—")
    assert (c[3], c[5], c[6]) == ("6.0000", "-50.4000", "-400.0000 / -10.8846")


def test_a_reference_row_totals_to_its_own_energy(store):
    table.drop_row(_row("uhf", "-", "-", e_method=-50.0069, e2_meh=None))
    c = _cells(_body(store)[0])
    assert (c[7], c[8], c[9]) == ("—", "-50.0069", "-6.9000 / -0.1878")


def test_render_refuses_broken_bookkeeping(store):
    bad = _row("casci", "free@pur", "(4,12)")
    bad["e_tot"] += 1e-3
    (store / "rows.json").write_text(json.dumps({table.row_key(bad): bad}))
    with pytest.raises(ValueError, match="violates E_tot = E_method"):
        table.render_table()


def test_conv_cell_forms():
    assert table.conv_of(3.2e-09, 11) == "3.2e-09/11"
    assert table.conv_of(3.5e-07, 200, converged=False) == "3.5e-07/200*"
    assert table.conv_of(1e-12, None) == "1.0e-12/—"
    assert table.conv_of(None, 7) == "—/7"


def test_fmt_range_collapses_to_one_value_when_the_set_agrees():
    assert table.fmt_range([72.0, 72.0000001]) == "72.0000"
    assert table.fmt_range([71.7847, 71.9654]) == "71.7847..71.9654"


def test_store_and_table_live_beside_the_legs():
    # the committed store and table are the family's records: anchored
    # on the package's parent (the legs' directory), wherever the
    # process runs from
    home = table.ROWS.parent
    assert table.TABLE.parent == home
    assert home.name == "u14e" and (home / "vv").is_dir()
    assert table.ROWS.name == "vv.json" and table.TABLE.name == "VV.md"
    assert table.ROWS.exists()
