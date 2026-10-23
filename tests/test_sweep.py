"""The window x treatment sweep: declaration, plan, quarantine, render.

Fast tier -- the runner's scheduling logic, no SCF."""

import math

import numpy as np
import pytest

from atomterm.sweep import (
    CellResult,
    SweepSpec,
    _wtag,
    plan_sweep,
    render_insert,
)


def _window(label, ncas, m_complete=True):
    from atomterm.windows import Window
    return Window(label, None, np.eye(ncas), m_complete)


def _spec(**kw):
    base = dict(label="test", ncore=1, na=2, nb=0,
                windows=(_window("w4", 4), _window("w5", 5)),
                treatments=("casci",))
    base.update(kw)
    return SweepSpec(**base)


def test_spectreatments_are_validated():
    with pytest.raises(AssertionError):
        _spec(treatments=("ccsd",))


def test_cell_keys_are_sanitized_and_unique():
    spec = _spec(windows=(_window("2f+2d", 4), _window("no 12 cols", 4)))
    tags = {_wtag(w) for w in spec.windows}
    assert "2f_2d" in tags and "no_12_cols" in tags
    keys = {spec.cell_key(t, "casci") for t in tags}
    assert len(keys) == 2
    assert all(k.startswith("sweep_test_") for k in keys)


def test_plan_is_the_manifest_before_anything_runs():
    text = plan_sweep(_spec(budgets={"sa": {"max_cycle_macro": 7}}))
    assert "SWEEP test" in text
    assert "dets = 6" in text            # C(4,2) = 6 for the 4-orbital window
    assert "m_complete = True" in text
    assert "max_cycle_macro" in text     # budgets are declared, not default
    assert "w4xcasci" in text            # the cell list is explicit


def test_render_insert_quarantines_errors():
    ok = CellResult(key="c_ok", status="ok", energy=-1.5, l2=2.0,
                    span_lz=1e-13, span_lperp=2e-13, dets=6)
    bad = CellResult(key="c_bad", status="error", error="boom: exploded")
    block = render_insert({"c_ok": vars(ok), "c_bad": vars(bad)})
    assert "<!-- sweep:begin -->" in block and "<!-- sweep:end -->" in block
    assert "c_ok | ok | -1.500000" in block
    assert "ERROR" in block and "boom" in block
    # the markers bracket exactly the generated content
    assert block.index("sweep:begin") < block.index("c_ok") < \
        block.index("sweep:end")
