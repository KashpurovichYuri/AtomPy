"""The trusted rungs (loaded, never re-run -- the cache is the
reproducibility unit), the reference energy, the committed anchors."""

import json
import sys
from pathlib import Path

import numpy as np  # noqa: E402
from pyscf import scf  # noqa: E402

from atomterm.cases import campaign_spec  # noqa: E402

PAYLOAD = Path("elements") / "U14e" / "perturbation.json"
RUNGS = Path("elements") / "U14e" / "rungs"


# the rungs -----------------------------------------------------------

def load_rung(chk: str) -> tuple:
    """One rung, LOADED (never re-run -- the cache is the
    reproducibility unit): the molecule from the chk, an SCF object
    rebuilt around it, the converged record spliced in."""
    from pyscf.lib import chkfile as _chk
    mol = _chk.load_mol(chk)
    rec = _chk.load(chk, "scf")
    if not isinstance(rec, dict) or "mo_coeff" not in rec:
        sys.exit(f"{chk}: no converged SCF record")
    uhf = np.asarray(rec["mo_coeff"]).ndim == 3
    mf = scf.UHF(mol) if uhf else scf.ROHF(mol)
    mf.__dict__.update(rec)
    mf.converged = True
    return mol, mf


def load_reference() -> tuple:
    """The trusted ROHF rung; the chk's identity scalars cross-checked
    against a freshly built declared molecule."""
    mol, mf = load_rung(str(RUNGS / "rohf.chk"))
    from atomterm.system import build_mol
    declared = build_mol(campaign_spec("U14e").atom_config())
    if (mol.nelectron != declared.nelectron or mol.spin != declared.spin
            or mol.nao != declared.nao):
        sys.exit(f"the ROHF rung's identity drifted from the "
                 f"declaration (nelectron {mol.nelectron}, spin "
                 f"{mol.spin}, nao {mol.nao} vs {declared.nelectron}/"
                 f"{declared.spin}/{declared.nao})")
    print(f"  ROHF (trusted rung, loaded): E = {mf.e_tot:.8f}, nao "
          f"{mol.nao}", flush=True)
    return mol, mf


_ROHF_E = None


def rohf_energy() -> float:
    """The reference energy -- the zero of every dE in the table."""
    global _ROHF_E
    if _ROHF_E is None:
        _, mf = load_rung(str(RUNGS / "rohf.chk"))
        _ROHF_E = float(mf.e_tot)
    return _ROHF_E


# the committed anchors ------------------------------------------------

def anchor(addr: str) -> dict:
    """The committed payload row an assertion reads -- fail fast
    unless it prices roots (the declared-roots discipline)."""
    payload = json.loads(PAYLOAD.read_text())
    for c in payload["cells"]:
        if c["address"] == addr:
            if c.get("disposition") != "rows" or not c.get("roots"):
                sys.exit(f"anchor {addr} carries no priced roots")
            return c
    sys.exit(f"anchor {addr} not in {PAYLOAD}")
