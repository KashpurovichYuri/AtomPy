"""The orbital-optimization engine: CASCI, CASSCF, and SA-CASSCF --
one deep driver (constrained_casscf) with declarative restriction
masks, plus the one-line rung adapters.

Conventions
-----------
- Active space: the active orbitals sit at ``[ncore : ncore + ncas]``
  of the FULL MO matrix -- CASCI/CASSCF slice it themselves; a
  truncated matrix fails ``nvir >= 0`` in ``check_sanity``.
- mcscf objects MUTATE: ``mc.kernel()`` canonicalizes
  ``mc.mo_coeff``, and ``get_h1eff()``/``get_h2eff()`` re-process
  whatever ``mc.mo_coeff`` currently holds -- pass ``mo`` explicitly
  and use a fresh object for integrals. This module's discipline: the
  SA roots pass takes ``roots.mo_coeff = mc.mo_coeff`` explicitly and
  reports the window EXACTLY as converged (``canonicalization``
  disabled).
- Checkpoints: every run persists energies + CI under
  ``atomterm/<name>/`` in the chkfile. Labels come from exact
  diagonalization of our H in the converged window -- PySCF's ci
  vectors straggle inside degenerate clusters and are never a label
  source.
- Mask vocabulary: every mask name starts with ``mask_`` and states
  the constraint (what is frozen / what never mixes / which rotations
  survive). The geometry masks are label-free and live here
  (``mask_frozen_window``, ``mask_no_core_window_mix``); the
  label-built ``windows.mask_radial_rotations`` completes the family.
  The serialized axis tokens of the cell addresses (``frozen``,
  ``virtmix``, ``radial``) are DATA identifiers, not names -- their
  mapping to these masks is stated once at the data layer (the cells
  docstring). The rung-purpose prose (core relaxation, virtual
  mixing) belongs to docstrings and README, anchored to these names.
"""

from __future__ import annotations

import types
from dataclasses import dataclass, field

import numpy as np
from pyscf import lib, mcscf

__all__ = [
    "CASResult",
    "run_casci",
    "constrained_casscf",
    "mask_frozen_window",
    "mask_no_core_window_mix",
    "run_casscf",
    "run_sa_casscf",
    "run_casscf_core_relax",
    "run_casscf_virtual_mix",
    "run_casscf_pyscf_start",
]


# ------------------------------------------------------------ the container

@dataclass
class CASResult:
    """Uniform result container for CASCI / CASSCF / SA-CASSCF.

    The stage's identity (tag, store name) flows SEPARATELY through
    stage_facts and the chkfile keys -- the container carries results
    only.

    Attributes
    ----------
    energies : np.ndarray
        Total energies of the roots (active-space ``e_cas`` + the
        common core part). For SA runs: the exact per-root CASCI pass
        in the converged window, never the ensemble value.
    ci : list
        CI vectors per root, PySCF order -- consumers indexing them
        against detbasis-order objects must map through
        ``cells.pyscf_ci_permutation`` (the ordering bridge).
    mo_coeff : np.ndarray or None
        The converged FULL MO matrix.
    mc : object or None
        The live pyscf object. Read ``converged`` from it, nothing
        else -- the mutation rule above.
    core_energy : float
        ``e_tot`` minus the (weight-averaged for SA) active energy:
        the constant every per-root total shares.
    """

    energies: np.ndarray  # total energies of the roots
    ci: list = field(default_factory=list)
    mo_coeff: np.ndarray | None = None
    mc: object | None = None  # live pyscf object
    core_energy: float = 0.0  # total minus active-space energy


# ------------------------------------------------- checkpointing + the kernel

def _checkpoint(chkfile: str, name: str, energies, ci) -> None:
    """Persist one run under ``atomterm/<name>/`` (energies + CI);
    ``chkfile=None`` disables. The object-array CI serialization does
    not round-trip for 2-D roots -- the loadable ``ci_stack`` form is
    the cell layer's own (it does not round-trip for 2-D roots)."""
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
    """Kernel + checkpoint + container. Everything read after
    ``mc.kernel()`` reads the CONVERGED orbitals -- the kernel
    canonicalizes ``mc.mo_coeff`` (the fact-6 mutation)."""
    mc.fcisolver.nroots = nroots
    mc.verbose = verbose
    mc.chkfile = chkfile
    mc.kernel()
    energies = _root_energies(mc)  # nroots=1 reproduces float(mc.e_tot)
    ci = [mc.ci] if nroots == 1 else list(mc.ci)
    _checkpoint(chkfile, name, energies, ci)
    return CASResult(energies, ci, mc.mo_coeff, mc, _core_energy(mc))


# --------------------------------------------------------------------- CASCI

def run_casci(mf, ncas, nelecas, mo_full, chkfile, nroots=1, verbose=0,
              name="casci", fci_opts: dict | None = None) -> CASResult:
    """CASCI in the given window (no orbital optimization).

    Parameters
    ----------
    mf, ncas, nelecas, mo_full, chkfile, nroots, verbose
        As :func:`constrained_casscf` (``mo_full`` mandatory here: the
        window is fully the caller's choice).
    name : str
        Chkfile store key -- keep it distinct per variation so
        differently-folded runs never collide.
    fci_opts : dict, optional
        Attributes set on the fcisolver before the kernel (e.g.
        ``{"pspace_size": ndet}`` forces pyscf's EXACT
        full-diagonalization branch -- Davidson straggles inside
        degenerate clusters: measured on U 2f2d/Ms=2, 48 roots: up to
        38 mEh within 13-17-fold multiplets while both exact spectra
        agree to 1e-7).
    """
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


def mask_frozen_window(self, nmo, ncore, ncas, frozen):
    """Allowed: ONLY inactive<->virtual rotations; every generator
    touching the active block is blocked -- the window is frozen
    outright (the relax rung; cells axis token ``frozen``).
    Signature: the ``uniq_var_indices`` hook, documented once in the
    seam comment above."""
    mask = np.zeros((nmo, nmo), dtype=bool)
    mask[ncore + ncas :, :ncore] = True
    return mask


def mask_no_core_window_mix(self, nmo, ncore, ncas, frozen):
    """Allowed: every rotation INVOLVING a virtual (virtual<->core and
    virtual<->active); blocked: active<->inactive alone -- no orbital
    ever trades places between the core and the active window (the
    virtmix rung; cells axis token ``virtmix``).
    Signature: the ``uniq_var_indices`` hook, documented once in the
    seam comment above."""
    mask = np.zeros((nmo, nmo), dtype=bool)
    mask[ncore + ncas :, : ncore + ncas] = True
    return mask


# ---------------------------------------------------------------- the driver

def constrained_casscf(mf, ncas, nelecas, chkfile, *, mo_full=None, mask=None,
                       weights=None, verbose=0, name="casscf",
                       mc_opts: dict | None = None,
                       fci_opts: dict | None = None) -> CASResult:
    """CASSCF with declaratively constrained orbital rotations -- the ONE
    deep driver every orbital-optimization rung goes through.

    Parameters
    ----------
    mf : the reference SCF object (integrals and the total convention).
    ncas, nelecas : active-space size and electron pair sector.
    chkfile : checkpoint location (or ``None`` to disable persisting).
    mo_full : np.ndarray, optional
        Starting MOs (e.g. raw shell columns via
        :func:`atomterm.windows.build_full_mo`); ``None`` = PySCF's
        default start (canonical MOs, canonical window) -- the
        start-variant axis.
    mask : callable, optional
        An ``uniq_var_indices``-shaped callable (see the masks above,
        or build one per window from AO labels); ``None`` = free
        optimization.
    weights : sequence of float, optional
        State-average over these weights. PySCF's SA kernel returns
        only the ENSEMBLE energy, so per-root energies come from a
        CASCI pass in the SA-converged window
        (``canonicalization = False``: the window reported EXACTLY as
        converged). The returned container pairs those exact per-root
        energies and core_energy with the SA kernel's own CI vectors
        (the variational object; cells persists it as ``ci_stack``)
        -- the roots pass's exact per-root CI is additionally stored
        under ``<name>_roots`` for inspection.
    name : str
        Chkfile store key under ``atomterm/`` -- keep it distinct per
        variation so differently-constrained runs never collide.
    mc_opts : dict, optional
        Attributes set on the mc object before the kernel (optimizer
        knobs: ``max_cycle_macro``, ``conv_tol``, ...).
    fci_opts : dict, optional
        Attributes set on the fcisolver before the kernel (see
        :func:`run_casci`).
    verbose : int
        Verbosity passed to the mc object.
    """
    mc = mcscf.CASSCF(mf, ncas, nelecas)
    if mo_full is not None:
        mc.mo_coeff = mo_full
    for k, v in (mc_opts or {}).items():
        setattr(mc, k, v)
    for k, v in (fci_opts or {}).items():
        setattr(mc.fcisolver, k, v)
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


# ----------------------------------------------------------------- adapters
# One line each: a named (start, mask, weights) recipe over the ONE
# driver -- a new rung is a new adapter or a direct constrained_casscf
# call, never a new kernel path.

def run_casscf(mf, ncas, nelecas, mo_full, chkfile, verbose=0) -> CASResult:
    """Free CASSCF from the given start.

    Parameters: as :func:`constrained_casscf` with ``mo_full``
    mandatory (raw-shell columns by convention), no mask.
    """
    return constrained_casscf(mf, ncas, nelecas, chkfile, mo_full=mo_full,
                              name="casscf", verbose=verbose)


def run_sa_casscf(mf, ncas, nelecas, mo_full, weights, chkfile, verbose=0,
                  mask=None) -> CASResult:
    """SA-CASSCF over the given state weights (primary sector). PySCF's SA
    kernel returns only the ensemble energy; per-root energies come from
    the CASCI pass in the converged window -- see
    :func:`constrained_casscf` (``weights``).

    Parameters: as :func:`constrained_casscf` plus ``weights`` (the
    state-average weights over the roots).
    """
    return constrained_casscf(mf, ncas, nelecas, chkfile, mo_full=mo_full,
                              weights=list(weights), mask=mask,
                              name="sa_casscf", verbose=verbose)


def run_casscf_core_relax(mf, ncas, nelecas, mo_full, chkfile, verbose=0,
                          mc_opts: dict = None) -> CASResult:
    """CASSCF with the active-window span frozen ("CASCI + core relaxation").

    Mask :func:`mask_frozen_window`: only inactive<->virtual
    rotations -- the core polarizes in response to the correlated CI
    density (real physics, typically ~mEh) but no orbital ever enters
    or leaves the window: state identity and window closure preserved
    by construction. Peculiarity: the relaxed core is slightly
    anisotropic; the analysis fold averages it, so the rung's
    ``max|dE|`` vs PySCF measures that anisotropy -- reported, not
    hidden.

    Parameters: as :func:`constrained_casscf` plus ``mc_opts``
    (per-run optimizer budgets).
    """
    return constrained_casscf(mf, ncas, nelecas, chkfile, mo_full=mo_full,
                              mask=mask_frozen_window, mc_opts=mc_opts,
                              name="casscf_relax", verbose=verbose)


def run_casscf_virtual_mix(mf, ncas, nelecas, mo_full, chkfile, verbose=0) -> CASResult:
    """CASSCF with core<->window exchange blocked ("virtmix" rung).

    Mask :func:`mask_no_core_window_mix`: active<->virtual and
    inactive<->virtual allowed, active<->inactive blocked -- the window
    borrows radial/polarization character from the virtuals (real
    dynamic correlation; appears as radial correlation in the natural
    orbitals) but the state's identity cannot change. relax < virtmix <
    free CASSCF brackets the contributions (README · Motivation).

    Parameters: as :func:`constrained_casscf`.
    """
    return constrained_casscf(mf, ncas, nelecas, chkfile, mo_full=mo_full,
                              mask=mask_no_core_window_mix,
                              name="casscf_virtmix", verbose=verbose)


def run_casscf_pyscf_start(mf, ncas, nelecas, chkfile, verbose=0) -> CASResult:
    """CASSCF from PySCF's own default guess (canonical ROHF MOs, active =
    the canonical window ``[ncore:ncore+ncas]``), NOT the raw shell columns.
    Isolates the start-dependence of the optimizer basin (compare with
    :func:`run_casscf` on the same system).

    Parameters: as :func:`constrained_casscf` with no ``mo_full`` --
    PySCF's default start IS the experiment.
    """
    return constrained_casscf(mf, ncas, nelecas, chkfile,
                              name="casscf_pyscf", verbose=verbose)


