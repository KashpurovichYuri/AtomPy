"""Regenerate tests/carbon_reference.py -- the FROZEN carbon ROHF.

One single-threaded draw (the cleanest, most canonical), stored as
plain source: Python float repr round-trips exactly, so loading the
module reproduces the SCF state bitwise. The unit tier runs NO SCF.
Run ONLY when the basis or the SCF defaults change -- the frozen draw
is the unit tier's determinism contract (certified by double-run
sha256), so a casual regen silently re-baselines every unit pin.

Usage:
    python scripts/generate_fixture_carbon.py
"""

from __future__ import annotations

import pathlib
import sys

import numpy as np
from pyscf import lib

lib.num_threads(1)

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from atomterm.system import AtomConfig, build_mol  # noqa: E402
from atomterm.reference import run_rohf  # noqa: E402
from atomterm.windows import select_shell_ao  # noqa: E402


def block(name: str, array: np.ndarray) -> str:
    """One numpy block as Python source (float repr -- exactly
    round-tripping)."""
    body = np.array2string(
        np.asarray(array), separator=", ", threshold=10 ** 6,
        max_line_width=79, formatter={"float_kind": lambda x: repr(float(x))},
    )
    return f"{name} = np.array({body}, dtype=float)"


def main() -> None:
    """One fresh single-threaded draw, written as
    tests/carbon_reference.py and certified bitwise on the way out."""
    chk = pathlib.Path("ab.tmp_freeze.chk")
    if chk.exists():
        chk.unlink()  # always a FRESH single-threaded draw
    cfg = AtomConfig(verbose=0)
    mol, mf = run_rohf(cfg, str(chk))
    assert mf.converged

    text = f'''"""The FROZEN carbon ROHF reference (C/cc-pVDZ, spin=2,
one single-threaded draw): the unit tier's deterministic substrate.
Floats are repr-round-trip exact -- loading this module reproduces the
SCF state bitwise and NO SCF runs in the unit tier. Regenerate (only
if the basis or defaults ever change) with
scripts/generate_fixture_carbon.py."""

import numpy as np

E_TOT = {float(mf.e_tot)!r}
{block("MO_ENERGY", mf.mo_energy)}
{block("MO_OCC", mf.mo_occ)}
{block("MO_COEFF", mf.mo_coeff)}
'''
    out = pathlib.Path("tests/carbon_reference.py")
    out.write_text(text, encoding="utf-8")
    chk.unlink(missing_ok=True)

    # certify the round trip bitwise before leaving
    sys.path.insert(0, "tests")
    import carbon_reference as ref  # noqa: E402
    assert ref.E_TOT == mf.e_tot
    assert ref.MO_COEFF.tobytes() == np.asarray(mf.mo_coeff).tobytes()
    assert ref.MO_ENERGY.tobytes() == np.asarray(mf.mo_energy).tobytes()
    assert ref.MO_OCC.tobytes() == np.asarray(mf.mo_occ).tobytes()
    print(f"froze {out} ({mol.nao_nr()} AOs) -- round trip bitwise")


if __name__ == "__main__":
    main()
