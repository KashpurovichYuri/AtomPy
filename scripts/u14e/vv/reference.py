"""The reference the routes run on.

The trusted rungs are LOADED, never re-run: a cached chkfile is the
reproducibility unit (a from-scratch SCF draw is chaotic).
:class:`Substrate` is the C1 reference -- the molecule and its SCF
record -- plus its D2h-tagged wrapper: the same AO quantities on a
molecule that carries the symmetry framework the pinned route's
solvers need.  No SCF ever runs in this module.

Layering
--------
``load_rung`` / ``load_reference`` / ``rohf_energy`` -- the trusted
rungs; ``anchor`` -- the committed payload rows;
``wrapper_certificate`` / ``wrapper_substrate`` -> :class:`Substrate`
-- the substrate.
"""

import json
from pathlib import Path

import numpy as np
from pyscf import gto, scf

from atomterm.cases import campaign_spec
from atomterm.system import build_mol

ROOT = Path(__file__).resolve().parents[3]


def rung_path(atom: str, name: str) -> Path:
    """The path of a trusted rung's chkfile (under ``elements/``)."""
    return ROOT / "elements" / atom / "rungs" / name


# ---------------------------------------------------------- the trusted rungs

def load_rung(chk) -> tuple:
    """One rung, LOADED (never re-run -- the cache is the
    reproducibility unit): the molecule from the chk, an SCF object
    rebuilt around it, the converged record spliced in.

    Parameters
    ----------
    chk : path of the chkfile.

    Returns
    -------
    (mol, mf) : the rung's molecule and its SCF object (ROHF, or UHF
        when the record's orbitals are spin-split).

    Raises
    ------
    RuntimeError : the chk carries no converged SCF record.
    """
    from pyscf.lib import chkfile as _chk
    mol = _chk.load_mol(str(chk))
    rec = _chk.load(chk, "scf")
    if not isinstance(rec, dict) or "mo_coeff" not in rec:
        raise RuntimeError(f"{chk}: no converged SCF record")
    uhf = np.asarray(rec["mo_coeff"]).ndim == 3
    mf = scf.UHF(mol) if uhf else scf.ROHF(mol)
    mf.__dict__.update(rec)
    mf.converged = True
    return mol, mf


def load_reference(atom: str = "U14e") -> tuple:
    """The trusted ROHF rung of a declared atom.

    The chk's identity scalars are cross-checked against a freshly
    built declared molecule (a drifted rung must never feed a leg).

    Returns
    -------
    (mol, mf) : the reference.

    Raises
    ------
    RuntimeError : the rung's identity drifted from the declaration.
    """
    mol, mf = load_rung(rung_path(atom, "rohf.chk"))
    declared = build_mol(campaign_spec(atom).atom_config())
    if (mol.nelectron != declared.nelectron or mol.spin != declared.spin
            or mol.nao != declared.nao):
        raise RuntimeError(
            f"the ROHF rung's identity drifted from the declaration "
            f"(nelectron {mol.nelectron}, spin {mol.spin}, nao "
            f"{mol.nao} vs {declared.nelectron}/{declared.spin}/"
            f"{declared.nao})")
    print(f"  ROHF (trusted rung, loaded): E = {mf.e_tot:.8f}, nao "
          f"{mol.nao}", flush=True)
    return mol, mf


_ROHF_E: dict[str, float] = {}


def rohf_energy(atom: str = "U14e") -> float:
    """The reference energy -- the zero of every dE in the table."""
    if atom not in _ROHF_E:
        _, mf = load_rung(rung_path(atom, "rohf.chk"))
        _ROHF_E[atom] = float(mf.e_tot)
    return _ROHF_E[atom]


# ------------------------------------------------------------- the anchors

def payload_path(atom: str = "U14e") -> Path:
    """The committed PT2 payload of a declared atom."""
    return ROOT / "elements" / atom / "perturbation.json"


def anchor(addr: str, atom: str = "U14e") -> dict:
    """The committed payload cell an assertion reads.

    Fails fast unless it prices roots (the declared-roots
    discipline): its ``roots`` carry the stored PT2 rows.

    Raises
    ------
    KeyError : the address is not in the payload.
    ValueError : the cell carries no priced roots.
    """
    payload = json.loads(payload_path(atom).read_text())
    for c in payload["cells"]:
        if c["address"] == addr:
            if c.get("disposition") != "rows" or not c.get("roots"):
                raise ValueError(f"anchor {addr} carries no priced roots")
            return c
    raise KeyError(f"anchor {addr} not in {payload_path(atom)}")


# ------------------------------------------------------ the D2h-tagged wrapper

def wrapper_certificate(mol_c1, mf_c1, mol_d2h, mf_d2h) -> dict:
    """Whether a D2h-tagged clone IS the C1 molecule: the same AO
    count and the same overlap and one-electron Hamiltonian exactly
    (AO integrals are gauge-independent; a difference means the clone
    is a different basis).

    Returns
    -------
    ``{"identical": bool, "nao": (n_c1, n_d2h), "dev_s": float,
    "dev_h": float}`` -- the deviations are ``max|A_c1 - A_d2h|``.
    """
    nao = (mol_c1.nao, mol_d2h.nao)
    dev_s = float(np.abs(mf_c1.get_ovlp() - mf_d2h.get_ovlp()).max()) \
        if nao[0] == nao[1] else float("inf")
    dev_h = float(np.abs(mf_c1.get_hcore() - mf_d2h.get_hcore()).max()) \
        if nao[0] == nao[1] else float("inf")
    return {"identical": nao[0] == nao[1] and dev_s == 0.0 and dev_h == 0.0,
            "nao": nao, "dev_s": dev_s, "dev_h": dev_h}


def wrapper_substrate(mol_c1, mf_c1, cfg) -> tuple:
    """The D2h-tagged wrapper of a C1 reference: the clone molecule
    (the atom, basis, spin and ECP ``cfg`` declares, tagged ``D2h``)
    around an SCF object carrying the C1 rung's record -- no SCF
    runs.  The CI consumes AO integrals (gauge-independent) plus
    given columns; the C1 rung supplies both, the tagged molecule
    supplies the symmetry framework (the block projectors, the
    ``wfnsym``-pinned solvers).

    Parameters
    ----------
    mol_c1, mf_c1 : the C1 reference.
    cfg : the atom's declaration (``symbol``, ``basis``, ``spin``,
        ``ecp`` -- ``system.AtomConfig`` or a case's ``atom_config()``).

    Returns
    -------
    (mol, mf) : the D2h-tagged molecule and the SCF object on it
        (``mf.mo_coeff`` IS the C1 record's matrix).

    Raises
    ------
    ValueError : the clone's AO quantities differ from the C1
        molecule's (:func:`wrapper_certificate`).
    """
    mol_d2h = gto.M(atom=f"{cfg.symbol} 0 0 0", basis=cfg.basis,
                    spin=cfg.spin, ecp=cfg.ecp, symmetry="D2h",
                    verbose=0)
    mf_d2h = scf.ROHF(mol_d2h)
    mf_d2h.__dict__.update({"mo_coeff": mf_c1.mo_coeff,
                            "mo_occ": mf_c1.mo_occ,
                            "e_tot": mf_c1.e_tot, "converged": True})
    cert = wrapper_certificate(mol_c1, mf_c1, mol_d2h, mf_d2h)
    if not cert["identical"]:
        raise ValueError(
            f"the D2h clone is not the C1 molecule (nao {cert['nao']}, "
            f"ovlp dev {cert['dev_s']:.2e}, hcore dev {cert['dev_h']:.2e}) "
            f"-- the declaration {cfg.symbol}/{cfg.basis!r} does not "
            f"describe the reference")
    return mol_d2h, mf_d2h


class Substrate:
    """The reference the routes run on.

    Attributes
    ----------
    mol, mf : the C1 molecule and its SCF record -- the free routes'
        substrate.
    spec : the campaign's case declaration (the windows read its
        shell layer); None when the substrate is not a declared case.

    ``d2h`` is the D2h-tagged wrapper ``(mol, mf)`` -- the pinned
        route's substrate: built once on first use from ``spec`` and
        certified there, or supplied ready-made (a reference that is
        already symmetry-tagged supplies itself).
    """

    def __init__(self, mol, mf, spec=None, d2h=None):
        self.mol, self.mf, self.spec = mol, mf, spec
        self._d2h = d2h

    @classmethod
    def load(cls, atom: str = "U14e") -> "Substrate":
        """The trusted reference of a declared atom."""
        mol, mf = load_reference(atom)
        return cls(mol, mf, spec=campaign_spec(atom))

    @property
    def d2h(self) -> tuple:
        if self._d2h is None:
            self._d2h = wrapper_substrate(self.mol, self.mf, self._config())
        return self._d2h

    def _config(self):
        if self.spec is None:
            raise ValueError(
                "the substrate carries no case declaration: supply the "
                "D2h wrapper (d2h=) or a spec to build it from")
        return self.spec.atom_config()
