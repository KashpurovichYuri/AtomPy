"""Calculation drivers: ROHF, CASCI, CASSCF, SA-CASSCF with .chkfile caching.

Every stage writes a PySCF checkpoint file (MOs, densities) so re-runs are
fast: if the ROHF chkfile exists it is loaded instead of recomputed. Analysis
artifacts (energies, CI vectors) are stored under the ``atomterm/`` key
prefix of the same chkfile.

Active-space convention: the active orbitals are the AO functions of the
target valence shell (e.g. 2p), placed at ``[ncore : ncore+ncas]`` of the full
MO matrix. Same-shell AOs are orthonormal and exactly rotation-closed, which
is what makes the L-quantum-number analysis exact.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from pyscf import lib, mcscf, scf

__all__ = [
    "run_rohf",
    "select_shell_ao",
    "purified_s_core",
    "complete_orthonormal",
    "build_full_mo",
    "run_casci",
    "run_casscf",
    "run_sa_casscf",
    "CASResult",
]


def run_rohf(cfg, chkfile: str):
    """ROHF for the configured atom, cached in ``chkfile``.

    If the chkfile exists, the converged MOs are loaded and no SCF is run.
    Returns ``(mol, mf)``.
    """
    from atomterm.system import build_mol

    mol = build_mol(cfg)
    try:
        _, rec = scf.chkfile.load_scf(chkfile)
        mf = scf.ROHF(mol)
        mf.__dict__.update(rec)
        mf.converged = True
        mf.verbose = cfg.verbose
        return mol, mf
    except (OSError, KeyError, ValueError):
        pass
    mf = scf.ROHF(mol)
    mf.chkfile = chkfile
    mf.verbose = cfg.verbose
    mf.conv_tol = 1e-12  # tight core density: keeps rotational symmetry clean
    mf.run()
    return mol, mf


def select_shell_ao(mol, pattern: str) -> list[int]:
    """AO indices of a named valence shell, e.g. pattern="2p" -> [3, 4, 5].

    Plain prefix matching on PySCF's AO labels (``search_ao_label`` patterns
    are regex-flavored and over-match: "2p*" also hits "2s").
    """
    labels = [line.split()[-1] for line in mol.ao_labels()]
    return [i for i, lbl in enumerate(labels) if lbl.startswith(pattern)]


def purified_s_core(mol, mf, n_core: int):
    """Canonical core MOs projected onto the s-AO subspace (l = 0 exactly).

    Why this is needed (causal chain):

    1. Numerically converged SCF orbitals of an open-shell atom are NOT exact
       symmetry eigenfunctions: the reference density (e.g. |px py| for C) is
       anisotropic, the converged Fock inherits that anisotropy, and the
       canonical "s" MOs carry a small p/d tail (measured on C/cc-pVDZ:
       1s 1.9e-7, 2s 6.3e-6 weight deficit). The tail is numerical noise,
       not physics -- the core of a central-field atom is spherical.
    2. That tail breaks two preconditions at once:
       - orthonormality: <MO_1 | 2p_AO> ~ 1e-3 != 0, so a W built from such
         columns is not orthonormal, and second quantization (which assumes
         an orthonormal one-particle basis) is ill-defined on it;
       - rotational closure: L maps a pure solid harmonic within its (n,l)
         shell only; L acting on the p-tail leaves the span, so [H, Lz]
         fails at the ~1e-4 level instead of machine precision.
    3. Projection onto the s-AO subspace is the orbital-space analog of the
       spherically averaged density used in the folded path: it takes the
       l = 0 component exactly. The energy shifts only at second order in
       the tail amplitude (~1e-10 Eh) -- noise repair, not a model change.
    4. A Loewdin orthonormalization inside the remaining pair then restores
       orthonormality while staying inside the l = 0 subspace (rotations
       within l = 0 preserve purity).

    Returns ``(C, contamination)``: full-length (nao, n_core) orthonormal
    columns and the fraction ``1 - <c_s|S|c>`` of each canonical MO that was
    dropped (the SCF symmetry-breaking diagnostic). Note the MO s-weight
    lives in the AO metric (c^dag S c = 1), not the Euclidean norm.
    """
    import re

    labels = [line.split()[-1] for line in mol.ao_labels()]
    s_idx = [i for i, lbl in enumerate(labels) if re.fullmatch(r"\d+s", lbl)]

    c_can = mf.mo_coeff[:, :n_core]
    c_s = c_can[s_idx, :]                                # drop the p/d tail
    # MO normalization lives in the AO metric (c^dag S c = 1), so the
    # s-subspace weight is c_s^dag S_ss c_s -- not the Euclidean norm.
    s_ov = mf.get_ovlp()[np.ix_(s_idx, s_idx)]
    weight = np.einsum("ij,ij->j", c_s, s_ov @ c_s)
    contamination = np.clip(1.0 - weight, 0.0, None)     # fp noise can dip < 0
    c_s = c_s / np.sqrt(weight)                          # renormalize the rest

    gram = c_s.conj().T @ s_ov @ c_s                     # Loewdin in S metric
    ev, u = np.linalg.eigh(gram)
    c_s = c_s @ (u / np.sqrt(ev)) @ u.T

    c = np.zeros((mol.nao_nr(), n_core), dtype=c_s.dtype)
    c[s_idx, :] = c_s
    return c, contamination


def complete_orthonormal(v: np.ndarray, s: np.ndarray) -> np.ndarray:
    """Append an orthonormal-in-metric-``s`` complement to the columns of ``v``."""
    ev, u = np.linalg.eigh(s)
    s_half = (u * np.sqrt(ev)) @ u.T
    s_ihalf = (u / np.sqrt(ev)) @ u.T
    vx = s_half @ v
    resid = np.eye(v.shape[0]) - vx @ vx.conj().T
    p, sv, _ = np.linalg.svd(resid)
    c = s_ihalf @ p[:, sv > 1e-8]
    return np.hstack([v, c])


def build_full_mo(mf, core_idx: list[int], act_idx: list[int]) -> np.ndarray:
    """Full orthonormal MO matrix with chosen core and active columns placed
    so that CASCI/CASSCF's window ``[ncore : ncore+ncas]`` selects the active
    space.

    Core columns are canonical SCF MOs; active columns are the RAW AO
    functions of the target shell (columns of the identity), which are exactly
    rotation-closed -- essential for the L-quantum-number analysis.
    """
    n = mf.mo_coeff.shape[0]
    core = mf.mo_coeff[:, core_idx]
    act = np.eye(n)[:, act_idx]
    v = np.hstack([core, act])
    if v.shape[1] == n:
        return v
    return complete_orthonormal(v, mf.get_ovlp())


@dataclass
class CASResult:
    """Uniform result container for CASCI / CASSCF / SA-CASSCF."""

    kind: str
    energies: np.ndarray  # total energies of the roots
    ci: list = field(default_factory=list)
    mo_coeff: np.ndarray = None
    mc: object = None  # live pyscf object
    core_energy: float = 0.0  # total minus active-space energy


def _persist(chkfile: str, name: str, energies, ci) -> None:
    if chkfile is None:
        return
    lib.chkfile.save(chkfile, f"atomterm/{name}/energies", np.asarray(energies))
    lib.chkfile.save(chkfile, f"atomterm/{name}/ci", np.asarray(ci, dtype=object))


def _core_energy(mc) -> float:
    """total(e_tot) minus active-space energy (weight-aware for SA runs)."""
    e_cas = np.atleast_1d(np.asarray(mc.e_cas)).ravel().astype(float)
    e_tot = float(np.atleast_1d(np.asarray(mc.e_tot)).ravel()[0])
    weights = getattr(mc.fcisolver, "weights", None)
    if weights is not None and len(weights) == len(e_cas):
        avg_active = float(np.dot(weights, e_cas))
    else:
        avg_active = float(e_cas[0])
    return e_tot - avg_active


def _root_energies(mc) -> np.ndarray:
    """Total energies of all roots (active-space e_cas + common core part)."""
    e_cas = np.atleast_1d(np.asarray(mc.e_cas)).ravel().astype(float)
    return e_cas + _core_energy(mc)


def _run(mc, nroots: int, chkfile: str, name: str, verbose: int) -> CASResult:
    mc.fcisolver.nroots = nroots
    mc.verbose = verbose
    mc.chkfile = chkfile
    mc.kernel()
    if nroots == 1:
        energies = np.atleast_1d(float(np.atleast_1d(np.asarray(mc.e_tot)).ravel()[0]))
    else:
        energies = _root_energies(mc)
    ci = [mc.ci] if nroots == 1 else list(mc.ci)
    _persist(chkfile, name, energies, ci)
    return CASResult(name, energies, ci, mc.mo_coeff, mc, _core_energy(mc))


def run_casci(mf, ncas, nelecas, mo_full, chkfile, nroots=1, verbose=0) -> CASResult:
    mc = mcscf.CASCI(mf, ncas, nelecas)
    mc.mo_coeff = mo_full
    return _run(mc, nroots, chkfile, "casci", verbose)


def run_casscf(mf, ncas, nelecas, mo_full, chkfile, verbose=0) -> CASResult:
    mc = mcscf.CASSCF(mf, ncas, nelecas)
    mc.mo_coeff = mo_full
    return _run(mc, 1, chkfile, "casscf", verbose)


def run_sa_casscf(mf, ncas, nelecas, mo_full, weights, chkfile, verbose=0) -> CASResult:
    mc = mcscf.CASSCF(mf, ncas, nelecas)
    mc.mo_coeff = mo_full
    mc.state_average_(list(weights))
    return _run(mc, len(weights), chkfile, "sa_casscf", verbose)
