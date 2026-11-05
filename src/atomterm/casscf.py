"""The orbital-optimization engine: CASCI, CASSCF, and SA-CASSCF --
one deep driver (constrained_casscf) with declarative restriction
masks, plus the classic pipeline's stage adapters (run_casscf,
run_sa_casscf -- the cases stage names; every other recipe is a direct
constrained_casscf call, as cells does).

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
  docstring). The treatment-purpose prose (core relaxation, virtual
  mixing) belongs to docstrings and README, anchored to these names.
"""

from __future__ import annotations

import types
from dataclasses import dataclass, field

import numpy as np
from pyscf import ao2mo, fci, lib, mcscf
from pyscf.mcscf import mc1step

from atomterm.irreps import branch_d2h, orb_irreps
from atomterm.operators import active_space_integrals
from atomterm.windows import build_full_mo, occupied_core_idx

__all__ = [
    "CASResult",
    "run_casci",
    "constrained_casscf",
    "mask_frozen_window",
    "mask_no_core_window_mix",
    "run_casscf",
    "run_sa_casscf",
    "run_irrep_casci",
    "run_irrep_sa",
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
        ``operators.pyscf_ci_permutation`` (the ordering bridge).
    mo_coeff : np.ndarray or None
        The converged FULL MO matrix.
    mc : object or None
        The live pyscf object. Read ``converged`` and the optimizer's
        ``cycles`` from it, nothing else -- the mutation rule above.
    core_energy : float
        ``e_tot`` minus the (weight-averaged for SA) active energy:
        the constant every per-root total shares.
    cycles : int or None
        The optimizer's macro-cycle count (cost pricing + restart
        audits); ``None`` for runs without an optimizer (CASCI).
    conv_dE : float or None
        The final macro-cycle |delta e_tot| -- the convergence
        precision; ``None`` without an optimizer or before the second
        cycle.
    """

    energies: np.ndarray  # total energies of the roots
    ci: list = field(default_factory=list)
    mo_coeff: np.ndarray | None = None
    mc: object | None = None  # live pyscf object
    core_energy: float = 0.0  # total minus active-space energy
    cycles: int | None = None  # optimizer macro-cycles (CI: None)
    conv_dE: float | None = None  # final |delta e_tot| (CI: None)


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


def _cycle_trace(mc) -> list:
    """Install the macro-cycle instrument on an optimizer and return
    its mutable cell -- PySCF never stores the count (the loop counter
    is kernel-local); the per-cycle ``callback`` IS the instrument,
    chained onto any user callback.  The cell reads ``[count, dE,
    e_prev]`` after ``kernel()``: the macro-cycle count, the final
    |delta e_tot| (the convergence precision; NaN before the second
    cycle), and the previous cycle's energy.
    """
    out = [0, float("nan"), None]
    prev = getattr(mc, "callback", None)

    def _tick(env):
        e = float(env.get("e_tot", float("nan")))
        if out[2] is not None:
            out[1] = abs(e - out[2])
        out[2] = e
        out[0] += 1
        if callable(prev):
            prev(env)

    mc.callback = _tick
    return out


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
    # the optimizer class only (mcscf.CASSCF is the factory;
    # mc1step.CASSCF the class; CASCI is not its subclass)
    trace = _cycle_trace(mc) if isinstance(mc, mc1step.CASSCF) else None
    mc.kernel()
    energies = _root_energies(mc)  # nroots=1 reproduces float(mc.e_tot)
    ci = [mc.ci] if nroots == 1 else list(mc.ci)
    _checkpoint(chkfile, name, energies, ci)
    return CASResult(
        energies, ci, mc.mo_coeff, mc, _core_energy(mc),
        cycles=trace[0] if trace is not None else None,
        conv_dE=(float(trace[1]) if trace is not None
                 and np.isfinite(trace[1]) else None))


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
        ``{"pspace_size": ndet}`` forces PySCF's exact
        full-diagonalization branch -- the cure when Davidson
        straggles inside degenerate multiplets).
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
    outright (core relaxation; cells axis token ``frozen``).
    Signature: the ``uniq_var_indices`` hook, documented once in the
    seam comment above."""
    mask = np.zeros((nmo, nmo), dtype=bool)
    mask[ncore + ncas :, :ncore] = True
    return mask


def mask_no_core_window_mix(self, nmo, ncore, ncas, frozen):
    """Allowed: every rotation INVOLVING a virtual (virtual<->core and
    virtual<->active); blocked: active<->inactive alone -- no orbital
    ever trades places between the core and the active window (virtual
    mixing; cells axis token ``virtmix``).
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
    deep driver every orbital-optimization treatment goes through.

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
        State-average over these weights. A ONE-root list IS the
        single root (the plain kernel -- PySCF's SA wrapper cannot
        einsum the scalar base return); longer lists run the SA
        kernel, which returns only the ENSEMBLE energy, so per-root
        energies come from a
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

    Raises
    ------
    ValueError
        ``weights=[]`` -- the empty list is a caller bug (the
        measured-block rule never yields it); refuse here, not
        inside PySCF's SA wrapper.
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
    if weights is not None and not len(weights):
        raise ValueError(
            "an empty weights list is a caller bug (the measured-block "
            "rule never yields []; pass weights=None for the plain "
            "single-root kernel)")
    if weights is None or len(weights) == 1:
        # a one-root SA IS the single root (PySCF's SA wrapper einsums
        # the energies -- with nroots=1 the base kernel returns a
        # scalar and the wrapper's 'i' subscript has no axis to walk)
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


# ------------------------------------------------------------ stage adapters
# One line each: a named (start, mask, weights) recipe over the ONE
# driver. The two survivors are the classic pipeline's stage names
# (cases); anything else calls constrained_casscf directly (cells) --
# never a new kernel path.

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
                              name="sa-casscf", verbose=verbose)




# ------------------------------------------------- the symmetry lane (pinned)
# PySCF-native evaluators over a term's D2h component set: ONE
# wfnsym-pinned solver per irrep of irreps.branch_d2h(L, parity), the
# multiplicity as its root count, on the substrate signature.  They
# are evaluators, not optimizer recipes: no mask or checkpoint
# vocabulary applies (nothing is cached; the cells axis arrives with
# the declaration wave), and the SA form IS the two-pass pattern --
# the ensemble mix kernel optimizes, run_irrep_casci re-solves the
# exact roots in the frozen converged window.

def _irrep_solvers(mol, orbsyms, L: int, parity: str, conv_tol: float):
    """``[(irrep, solver, multiplicity)]`` -- one pinned solver per
    component block of the term."""
    solvers = []
    for irrep, n in branch_d2h(L, parity).items():
        s = fci.direct_spin1_symm.FCI(mol)
        s.orbsym = np.asarray(orbsyms)
        s.wfnsym = irrep
        s.nroots = n
        s.conv_tol = conv_tol
        solvers.append((irrep, s, n))
    return solvers


def run_irrep_casci(mol, mf, c_act, ncore: int, na: int, nb: int,
                    L: int, parity: str, conv_tol: float = 1e-12,
                    c_core=None) -> dict:
    """The pinned route at the CI level: one block solve per component.

    Every irrep of the term's branching gets a ``wfnsym``-pinned
    solver with ``nroots`` = multiplicity, on the SAME folded h1 as
    our determinant-basis matrices (:func:`constrained_casscf`'s
    convention twin) and PySCF's own window eri feed.

    Parameters
    ----------
    mol, mf, c_act, ncore, na, nb : the substrate signature (``mol``
        must be a D2h molecule; ``c_act`` the window columns).
    L, parity : the target TERM -- ``parity`` is the configuration's
        total inversion parity (p^2 -> g, s p -> u), not the shell's.
    conv_tol : the pinned solvers' CI convergence tolerance.
    c_core : explicit fold columns (the ladder's non-prefix folds --
        ``ncore`` is then the fold's column count; the default None
        folds the ncore prefix).

    Returns
    -------
    ``{irrep: (e_cas roots, ci vectors in PySCF order)}`` -- the
    lowest ``multiplicity`` roots per block.  That the lowest roots
    ARE the target term's components is an assumption the caller
    verifies with the exact meters (analysis bridges the vectors).

    Raises
    ------
    ValueError : via :func:`atomterm.irreps.orb_irreps`, if the
        molecule is not D2h.
    """
    h1, _ = active_space_integrals(mol, mf, c_act, ncore, c_core=c_core)
    eri = ao2mo.kernel(mol, c_act)
    orbsyms = orb_irreps(mol, c_act)
    out = {}
    for irrep, solver, _ in _irrep_solvers(mol, orbsyms, L, parity,
                                           conv_tol):
        e, ci = solver.kernel(h1, eri, c_act.shape[1], (na, nb))
        out[irrep] = (np.atleast_1d(np.asarray(e, dtype=float)),
                      np.asarray(ci))
    return out


def run_irrep_sa(mol, mf, c_act, ncore: int, na: int, nb: int,
                 L: int, parity: str, conv_tol: float = 1e-12,
                 c_core=None) -> dict:
    """The pinned route at the SA level: the component ensemble, optimized.

    ``mcscf.state_average_mix_`` over the pinned solvers with uniform
    component weights (the free-atom ensemble over the ``2L + 1``
    components), then the exact per-component roots re-solved in the
    frozen converged window (:func:`run_irrep_casci` -- PySCF's mix
    kernel returns the ensemble only; the two-pass pattern).

    Parameters
    ----------
    As :func:`run_irrep_casci` (the optimizer convergence is PySCF's
    default; the readout solvers take ``conv_tol``), including the
    ``c_core`` fold passthrough.

    Returns
    -------
    ``{"e_tot": ensemble total, "mo": the converged FULL MO matrix,
    "converged": bool, "cycles": the optimizer's macro-cycle count,
    "conv_dE": the final macro-cycle |delta e_tot|,
    "roots": the run_irrep_casci readout}``.

    Notes
    -----
    The optimized window stays irrep-aligned (the symmetry-adapted
    optimizer rotates only within irrep blocks), so the readout's
    block structure survives optimization.  Never bit-compare ``mo``
    against the free route -- the two optimizer paths land in different
    degenerate orbital gauges; energies and closure certificates are
    the comparison surface.
    """
    ncas = c_act.shape[1]
    if c_core is None:
        mo_full = build_full_mo(mf, occupied_core_idx(mf, ncore), c_act)
    else:
        nf = c_core.shape[1]
        mo_full = build_full_mo(mf, list(range(nf)), c_act, c_core=c_core)
        ncore = nf
    active = slice(ncore, ncore + ncas)
    mc = mcscf.CASSCF(mf, ncas, (na, nb))
    mc.verbose = 0
    solvers = _irrep_solvers(mol, orb_irreps(mol, mo_full[:, active]),
                             L, parity, conv_tol)
    mcscf.state_average_mix_(
        mc, [s for _, s, _ in solvers], [1.0 / (2 * L + 1)] * (2 * L + 1))
    trace = _cycle_trace(mc)
    mc.kernel(mo_full)
    roots = run_irrep_casci(mol, mf, mc.mo_coeff[:, active], ncore,
                            na, nb, L, parity, conv_tol, c_core=c_core)
    return {"e_tot": float(mc.e_tot), "mo": mc.mo_coeff,
            "converged": bool(mc.converged), "cycles": trace[0],
            "conv_dE": (float(trace[1]) if np.isfinite(trace[1])
                        else None),
            "roots": roots}
