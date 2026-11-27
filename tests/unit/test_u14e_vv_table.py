"""The V&V table module (scripts/u14e/vv/table.py) on a temporary store
(the committed one is never touched): the store accepts only the
declared keys and only records in the schema, the table renders one
block per method class with a rule row between blocks and the plain
ladder before the @r1 rows, a row's total is COMPUTED (never
supplied), the withdrawn cells travel together, a withdrawn row
renders dashes in the withdrawn cells only, and the cells' forms are
the table's conventions."""

import json

import pytest

from vv import table
from vv.table import Route, Row

pytestmark = pytest.mark.fast


@pytest.fixture
def store(tmp_path, monkeypatch):
    """The row store and the table redirected to tmp_path; the ROHF
    energy (the zero of every dE) fixed."""
    monkeypatch.setattr(table, "ROWS", tmp_path / "rows.json")
    monkeypatch.setattr(table, "TABLE", tmp_path / "table.md")
    monkeypatch.setattr(table, "rohf_energy", lambda: -50.0)
    return tmp_path


def _row(method, route, cas, e_method=-50.4, e2=-0.3, l2="72.0000"):
    return Row(method, route, cas, e_method, e2, l2, "6.0000")


def _body(tmp_path):
    """The data rows and the rule rows (the title, the blank line, the
    header and its rule are skipped)."""
    return (tmp_path / "table.md").read_text().splitlines()[4:]


def _cells(line):
    return [c.strip() for c in line.split("|")[1:-1]]


# ------------------------------------------------------------- the vocabulary

def test_the_keys_are_the_keep_set_in_display_order():
    keys = table.table_keys()
    assert keys[:2] == ("rohf|-|-", "uhf|-|-")
    assert len(keys) == len(set(keys)) == 38
    assert "sa-casscf|free@pur|(4,12)" in keys
    assert "sa-casscf|pinned@pur|(4,12)@r1" in keys
    assert "casci|pinned@pur|(4,12)" not in keys       # not a kept route
    assert keys.index("sa-casscf|free@pur|(6,13)@r1") < \
        keys.index("sa-casscf|pinned@pur|(4,12)")       # free rows first
    assert [r.value for r in Route] == ["-", "free@pur", "pinned@pur"]


# ------------------------------------------------------------------ the row

def test_the_total_is_computed_never_supplied():
    row = _row("casci", Route.FREE, "(4,12)", e_method=-50.25, e2=-0.44)
    assert row.e_total == pytest.approx(-50.69)
    assert row.key == "casci|free@pur|(4,12)"
    assert "e_tot" not in row.to_json() and set(row.to_json()) == {
        "method", "route", "cas", "e_method", "e2", "l2", "s2", "conv"}
    assert Row.from_json(row.to_json()) == row


def test_a_reference_row_totals_to_its_own_energy_and_never_prices():
    row = Row("uhf", Route.REFERENCE, "-", -50.0069, None, "-", "6.2289 (det)")
    assert row.e_total == -50.0069 and row.reference and not row.withdrawn
    with pytest.raises(ValueError, match="reference row .* carries an E"):
        Row("rohf", Route.REFERENCE, "-", -50.0, -0.3, "-", "-")


def test_the_withdrawn_cells_travel_together():
    row = Row("sa-casscf", Route.PINNED, "(6,13)", -50.45, None, None, "12.0000",
              "7.6e-10/13")
    assert row.withdrawn and row.e_total is None
    for e2, l2 in ((None, "72.0000"), (-0.3, None)):
        with pytest.raises(ValueError, match="withdrawn together"):
            Row("sa-casscf", Route.PINNED, "(6,13)", -50.45, e2, l2, "12.0000")


def test_a_record_outside_the_schema_is_refused():
    good = _row("casci", Route.FREE, "(4,12)").to_json()
    with pytest.raises(ValueError, match="carries fields \\['e_tot'\\]"):
        Row.from_json({**good, "e_tot": -50.7})
    with pytest.raises(ValueError, match="carries fields \\['conv'\\]"):
        Row.from_json({k: v for k, v in good.items() if k != "conv"})
    with pytest.raises(ValueError):
        Row.from_json({**good, "route": "free"})       # the old spelling


# ----------------------------------------------------------- the store + table

def test_drop_row_refuses_an_undeclared_key(store):
    for row in (_row("casci", Route.FREE, "(4,99)"),            # off the ladder
                _row("casci", Route.PINNED, "(4,12)"),         # not a kept route
                _row("sa-casscf", Route.PINNED, "(4,18)")):
        with pytest.raises(ValueError, match="outside the declared table"):
            table.drop_row(row)
    assert not (store / "rows.json").exists()             # nothing was written


def test_render_refuses_an_undeclared_key_or_a_mismatched_record(store):
    rec = _row("casci", Route.FREE, "(4,12)").to_json()
    (store / "rows.json").write_text(json.dumps({"casci|free|(4,12)": rec}))
    with pytest.raises(ValueError, match=r"casci\|free\|\(4,12\)"):
        table.render_table()
    (store / "rows.json").write_text(json.dumps({"casci|free@pur|(4,15)": rec}))
    with pytest.raises(ValueError, match="disagrees with"):
        table.render_table()


def test_one_block_per_method_class_with_a_rule_between(store):
    for m in ("sa-casscf", "casscf", "casci"):             # dropped out of order
        table.drop_row(_row(m, Route.FREE, "(4,12)"))
    table.drop_row(_row("sa-casscf", Route.PINNED, "(4,12)"))
    table.drop_row(Row("rohf", Route.REFERENCE, "-", -50.0, None, "-", "-"))
    body = _body(store)
    assert [_cells(r)[0] for r in body if "===" not in r] == [
        "rohf", "casci (free@pur)", "casscf (free@pur)",
        "sa-casscf (free@pur)", "sa-casscf (pinned@pur)"]
    assert [i for i, r in enumerate(body) if "===" in r] == [1, 3, 5]


def test_within_a_block_the_plain_ladder_precedes_the_r1_rows(store):
    for cas in ("(6,13)@r1", "(4,16)", "(4,12)@r1", "(4,15)"):
        table.drop_row(_row("casci", Route.FREE, cas))
    assert [_cells(r)[1] for r in _body(store)] == [
        "(4,15)", "(4,16)", "(4,12)@r1", "(6,13)@r1"]   # size-ordered
    # the store file keeps the table's order, whatever the drop order
    assert [r.key for r in table.load_rows().values()] == [
        f"casci|free@pur|{c}" for c in ("(4,15)", "(4,16)", "(4,12)@r1", "(6,13)@r1")]


def test_a_withdrawn_row_renders_dashes_in_the_withdrawn_cells_only(store):
    table.drop_row(Row("sa-casscf", Route.PINNED, "(6,13)", -50.45, None, None,
                       "12.0000", "7.6e-10/13"))
    c = _cells(_body(store)[0])
    assert (c[2], c[7], c[8], c[9]) == ("—", "—", "—", "—")
    assert (c[3], c[4], c[5], c[6]) == ("12.0000", "7.6e-10/13", "-50.4500",
                                        "-450.0000 / -12.2451")


def test_a_reference_row_renders_its_own_total(store):
    table.drop_row(Row("uhf", Route.REFERENCE, "-", -50.0069, None, "-",
                       "6.2289 (det)"))
    c = _cells(_body(store)[0])
    assert (c[0], c[1], c[2], c[3], c[7], c[8], c[9]) == (
        "uhf", "-", "-", "6.2289 (det)", "—", "-50.0069", "-6.9000 / -0.1878")


def test_an_empty_store_renders_the_header_only(store):
    table.render_table()
    assert len((store / "table.md").read_text().splitlines()) == 4


def test_drop_row_merges_into_the_existing_records(store):
    table.drop_row(_row("casci", Route.FREE, "(4,15)"))
    table.drop_row(_row("casscf", Route.FREE, "(4,12)", e_method=-50.43))
    table.drop_row(_row("casci", Route.FREE, "(4,15)", e_method=-49.45))
    rows = table.load_rows()
    assert list(rows) == ["casci|free@pur|(4,15)", "casscf|free@pur|(4,12)"]
    assert rows["casci|free@pur|(4,15)"].e_method == -49.45   # replaced


def test_store_and_table_live_beside_the_legs():
    # the committed store and table are the family's records: anchored
    # on the package's parent (the legs' directory), wherever the
    # process runs from
    home = table.ROWS.parent
    assert table.TABLE.parent == home
    assert home.name == "u14e" and (home / "vv").is_dir()
    assert table.ROWS.name == "vv.json" and table.TABLE.name == "VV.md"
    assert table.ROWS.exists()


# ------------------------------------------------------------------ the cells

def test_conv_cell_forms():
    assert table.conv_of(3.2e-09, 11) == "3.2e-09/11"
    assert table.conv_of(3.5e-07, 200, converged=False) == "3.5e-07/200*"
    assert table.conv_of(1e-12, None) == "1.0e-12/—"
    assert table.conv_of(None, 7) == "—/7"


def test_fmt_range_collapses_to_one_value_when_the_set_agrees():
    assert table.fmt_range([72.0, 72.0000001]) == "72.0000"
    assert table.fmt_range([71.7847, 71.9654]) == "71.7847..71.9654"


def test_de_fmt_shows_both_units():
    assert table.de_fmt(-0.0069) == "-6.9000 / -0.1878"
