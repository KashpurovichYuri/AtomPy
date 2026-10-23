"""The orbital-optimization engine: CASCI, CASSCF, and SA-CASSCF --
one deep driver (constrained_casscf) with declarative restriction masks,
plus the one-line rung adapters.

Active-space convention: the active orbitals sit at
``[ncore : ncore + ncas]`` of the full MO matrix. Every stage writes a
checkpoint (energies and CI vectors under the ``atomterm/`` key prefix).
Labels come from exact diagonalization in the converged window -- the
PySCF ci vectors straggle inside degenerate clusters and are never a
label source.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from pyscf import lib, mcscf

__all__ = [
    "CASResult",
    "run_casci",
    "constrained_casscf",
    "core_virtual_only_mask",
    "no_core_window_exchange_mask",
    "run_casscf",
    "run_sa_casscf",
    "run_casscf_core_relax",
    "run_casscf_virtual_mix",
    "run_casscf_pyscf_start",
]


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
    energies = _root_energies(mc)  # nroots=1 reproduces float(mc.e_tot)
    ci = [mc.ci] if nroots == 1 else list(mc.ci)
    _persist(chkfile, name, energies, ci)
    return CASResult(name, energies, ci, mc.mo_coeff, mc, _core_energy(mc))


def run_casci(mf, ncas, nelecas, mo_full, chkfile, nroots=1, verbose=0,
              name="casci", fci_opts: dict = None) -> CASResult:
    """CASCI in the given window. ``name`` is the chkfile store key -- keep
    it distinct per variation so differently-folded runs never collide.
    ``fci_opts``: attributes set on the fcisolver before the kernel (e.g.
    ``{"pspace_size": ndet}`` forces pyscf's EXACT full-diagonalization
    branch -- Davidson straggles inside degenerate clusters: measured on
    U 2f2d/Ms=2, 48 roots: up to 38 mEh within 13-17-fold multiplets
    while both exact spectra agree to 1e-7)."""
    mc = mcscf.CASCI(mf, ncas, nelecas)
    mc.mo_coeff = mo_full
    for k, v in (fci_opts or {}).items():
        setattr(mc.fcisolver, k, v)
    return _run(mc, nroots, chkfile, name, verbose)


# ---------------------------------------------------------------------------
# Orbital-rotation masks (seam of :func:`constrained_casscf`).
#
# PySCF's CASSCF optimizer enumerates its rotation variables through the
# ``uniq_var_indices(nmo, ncore, ncas, frozen)`` hook, which returns a
# boolean (nmo, nmo) pair mask: the generator entry (i, j) is optimized
# iff mask[i, j]. Rebinding that method constrains the optimization WITHOUT
# touching the optimizer -- each mask below is one declarative constraint.
# ---------------------------------------------------------------------------


def core_virtual_only_mask(self, nmo, ncore, ncas, frozen):
    """ONLY inactive<->virtual rotations allowed (the relax rung)."""
    mask = np.zeros((nmo, nmo), dtype=bool)
    mask[ncore + ncas :, :ncore] = True
    return mask


def no_core_window_exchange_mask(self, nmo, ncore, ncas, frozen):
    """Everything EXCEPT active<->inactive rotations: no orbital ever
    trades places between the core and the active window (the virtmix
    rung)."""
    mask = np.zeros((nmo, nmo), dtype=bool)
    mask[ncore + ncas :, : ncore + ncas] = True
    return mask


def constrained_casscf(mf, ncas, nelecas, chkfile, *, mo_full=None, mask=None,
                       weights=None, verbose=0, name="casscf",
                       mc_opts: dict = None) -> CASResult:
    """CASSCF with declaratively constrained orbital rotations -- the one
    deep driver every orbital-optimization rung goes through.

    mo_full : starting MOs (e.g. raw shell columns via :func:`build_full_mo`);
        ``None`` = PySCF's default start (canonical MOs, canonical window) --
        the start-variant axis.
    mask : an ``uniq_var_indices``-shaped callable (see the masks above, or
        build one per window from AO labels); ``None`` = free optimization.
    weights : state-average over these weights. PySCF's SA kernel returns
        only the ENSEMBLE energy, so per-root energies come from a CASCI
        pass in the SA-converged window.
    name : chkfile store key under ``atomterm/`` -- keep it distinct per
        variation so differently-constrained runs never collide.
    mc_opts : attributes set on the mc object before the kernel (optimizer
        knobs: ``max_cycle_macro``, ``conv_tol``, ...).
    """
    import types

    mc = mcscf.CASSCF(mf, ncas, nelecas)
    if mo_full is not None:
        mc.mo_coeff = mo_full
    for k, v in (mc_opts or {}).items():
        setattr(mc, k, v)
    if mask is not None:
        mc.uniq_var_indices = types.MethodType(mask, mc)
    if weights is None:
        return _run(mc, 1, chkfile, name, verbose)
    mc.state_average_(list(weights))
    res = _run(mc, len(weights), chkfile, name, verbose)
    roots = mcscf.CASCI(mf, ncas, nelecas)
    roots.mo_coeff = mc.mo_coeff
    roots.canonicalization = False  # report the SA window exactly as converged
    rr = _run(roots, len(weights), chkfile, f"{name}_roots", verbose)
    res.energies = rr.energies
    res.core_energy = rr.core_energy
    return res


def run_casscf(mf, ncas, nelecas, mo_full, chkfile, verbose=0) -> CASResult:
    """Free CASSCF from the given start (raw-shell columns by convention)."""
    return constrained_casscf(mf, ncas, nelecas, chkfile, mo_full=mo_full,
                              name="casscf", verbose=verbose)


def run_sa_casscf(mf, ncas, nelecas, mo_full, weights, chkfile, verbose=0,
                  mask=None) -> CASResult:
    """SA-CASSCF over the given state weights (primary sector). PySCF's SA
    kernel returns only the ensemble energy; per-root energies come from
    the CASCI pass in the converged window -- see :func:`constrained_casscf`
    (``weights``)."""
    return constrained_casscf(mf, ncas, nelecas, chkfile, mo_full=mo_full,
                              weights=list(weights), mask=mask,
                              name="sa_casscf", verbose=verbose)


def run_casscf_core_relax(mf, ncas, nelecas, mo_full, chkfile, verbose=0) -> CASResult:
    """CASSCF with the active-window span frozen ("CASCI + core relaxation").

    Mask :func:`core_virtual_only_mask`: only inactive<->virtual
    rotations -- the core polarizes in response to the correlated CI
    density (real physics, typically ~mEh) but no orbital ever enters
    or leaves the window: state identity and window closure preserved
    by construction. Peculiarity: the relaxed core is slightly
    anisotropic; the analysis fold averages it, so the rung's
    ``max|dE|`` vs PySCF measures that anisotropy -- reported, not
    hidden.
    """
    return constrained_casscf(mf, ncas, nelecas, chkfile, mo_full=mo_full,
                              mask=core_virtual_only_mask,
                              name="casscf_relax", verbose=verbose)


def run_casscf_virtual_mix(mf, ncas, nelecas, mo_full, chkfile, verbose=0) -> CASResult:
    """CASSCF with core<->window exchange blocked ("virtmix" rung).

    Mask :func:`no_core_window_exchange_mask`: active<->virtual and
    inactive<->virtual allowed, active<->inactive blocked -- the window
    borrows radial/polarization character from the virtuals (real
    dynamic correlation; appears as radial correlation in the natural
    orbitals) but the state's identity cannot change. relax < virtmix <
    free CASSCF brackets the contributions (README §3.3).
    """
    return constrained_casscf(mf, ncas, nelecas, chkfile, mo_full=mo_full,
                              mask=no_core_window_exchange_mask,
                              name="casscf_virtmix", verbose=verbose)


def run_casscf_pyscf_start(mf, ncas, nelecas, chkfile, verbose=0) -> CASResult:
    """CASSCF from PySCF's own default guess (canonical ROHF MOs, active =
    the canonical window ``[ncore:ncore+ncas]``), NOT the raw shell columns.
    Isolates the start-dependence of the optimizer basin (compare with
    :func:`run_casscf` on the same system)."""
    return constrained_casscf(mf, ncas, nelecas, chkfile,
                              name="casscf_pyscf", verbose=verbose)


