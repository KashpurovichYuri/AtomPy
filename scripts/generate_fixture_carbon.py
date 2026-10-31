"""Regenerate the frozen carbon references -- the unit tier's
deterministic substrates.

One single-threaded draw (the cleanest, most canonical), stored as
plain source: Python float repr round-trips exactly, so loading the
module reproduces the SCF state bitwise. The unit tier runs NO SCF.
Run ONLY when the basis or the SCF defaults change -- a frozen draw
is the unit tier's determinism contract (certified by double-run
sha256), so a casual regen silently re-baselines every unit pin.

Usage:
    python scripts/generate_fixture_carbon.py          # the C1 ROHF
    python scripts/generate_fixture_carbon.py d2h      # the D2h twin
"""

from __future__ import annotations

import pathlib
import sys

import numpy as np
from pyscf import gto, lib, scf

lib.num_threads(1)

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from atomterm.system import AtomConfig  # noqa: E402
from atomterm.reference import run_rohf  # noqa: E402


def block(name: str, array: np.ndarray) -> str:
    """One numpy block as Python source (float repr -- exactly
    round-tripping)."""
    body = np.array2string(
        np.asarray(array), separator=", ", threshold=10 ** 6,
        max_line_width=79, formatter={"float_kind": lambda x: repr(float(x))},
    )
    return f"{name} = np.array({body}, dtype=float)"


def _certify(path: pathlib.Path, mf) -> None:
    """The fresh file must load back bitwise before we leave."""
    sys.path.insert(0, "tests")
    stem = path.stem
    ref = __import__(stem)
    assert ref.E_TOT == mf.e_tot
    assert ref.MO_COEFF.tobytes() == np.asarray(mf.mo_coeff).tobytes()
    assert ref.MO_ENERGY.tobytes() == np.asarray(mf.mo_energy).tobytes()
    assert ref.MO_OCC.tobytes() == np.asarray(mf.mo_occ).tobytes()
    print(f"froze {path} ({mf.mol.nao_nr()} AOs) -- round trip bitwise")


def main() -> None:
    """One fresh single-threaded C1 draw, written as
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
    _certify(out, mf)


def d2h() -> None:
    """One fresh single-threaded D2h draw, written as
    tests/carbon_reference_d2h.py (the symmetry-lane twin of the C1
    fixture -- same atom, symmetry-adapted; the reference is the
    deterministic D2h-symmetric solution, whatever its sphere
    breaking: the m-complete window re-symmetrizes at the CI)."""
    cfg = AtomConfig(verbose=0)
    mol = gto.M(atom=f"{cfg.symbol} 0 0 0", basis=cfg.basis, spin=cfg.spin,
                symmetry="D2h", verbose=0)
    mf = scf.ROHF(mol)
    mf.verbose = 0
    mf.kernel()
    assert mf.converged

    text = f'''"""The FROZEN carbon D2h twin (C/cc-pVDZ, spin=2,
symmetry-adapted single-threaded draw): the symmetry lane's
deterministic substrate. Floats are repr-round-trip exact -- loading
this module reproduces the SCF state bitwise and NO SCF runs in the
unit tier. Regenerate (only if the basis or defaults ever change)
with scripts/generate_fixture_carbon.py d2h."""

import numpy as np

E_TOT = {float(mf.e_tot)!r}
{block("MO_ENERGY", mf.mo_energy)}
{block("MO_OCC", mf.mo_occ)}
{block("MO_COEFF", mf.mo_coeff)}
'''
    out = pathlib.Path("tests/carbon_reference_d2h.py")
    out.write_text(text, encoding="utf-8")
    _certify(out, mf)


if __name__ == "__main__":
    {"d2h": d2h}.get(sys.argv[1] if len(sys.argv) > 1 else "", main)()
