"""Perturbation theory: the PT2 analyzation over certified CAS solutions.

PURPOSE
-------
Dynamic-correlation rows on top of any classified solution -- CASCI,
CASSCF, SA-CASSCF, any window.  The PT2 algebra is PySCF's own
``pyscf.mrpt`` (strongly contracted NEVPT2); this module owns the
BOUNDARY: it rebuilds the multi-root CASCI readout on the supplied
orbitals, refuses any drift from the stored solution (the bit-exact
assertion), and renders per-root, per-Koopmans-class rows.  The
reading discipline -- the gauge ledger, the two probes, the per-L
expectations -- lives in README ("Analysis"); the wave contract in
``wayfinder/specs/pt2-perturbation.md``.

Conventions
-----------
*   The substrate is ``(mol, mf, mo_full)`` + the CAS shape
    ``(ncas, nelecas)``.  ``mo_full``'s column layout IS the
    partition: CASCI freezes its own ``ncore`` prefix (derived from
    the electron counts -- the ECP-honest convention), keeps
    ``[ncore : ncore + ncas]`` active, the rest virtual.  The readout
    here is :func:`atomterm.casscf.run_casci`'s twin without the
    checkpoint: the same kernel call the stored cells were read out
    on, so a stored cell's energies recompute to machine zero.
*   Roots are indices into the multi-root readout (PySCF order),
    declared explicitly -- never re-derived, never inferred.
*   The stored-roots assertion FAILS FAST: a drifted readout refuses
    the whole call (``ValueError``) -- PT2 never runs on a solution
    the store does not own.  Vectors compare gauge-aware: an exactly
    degenerate block may lawfully re-solve in a rotated basis; the
    ROOT SET's projector must still match.
*   Per-root flow: the ``mrpt.NEVPT(mc, root=k)`` canonicalization
    prelude (``cas_natorb``; the kernel's own recipe), then the eight
    Koopmans-class solves through mrpt's public class functions; the
    class sum is the row's E(2), pinned against ``NEVPT(...).kernel()``
    in the tests (the plumbing-vs-kernel twin).  ``density_fit`` is
    pinned False: plain CASCI feeds the exact 4-index path.
*   The variant axis rides the row schema (``PT2Row.variant``);
    further PT2 families (CASPT2 on an env upgrade) join as new
    constructors at the same seam -- the registry is deferred to the
    second consumer (a seam with one consumer stays hypothetical).

Layering: pyscf only (mcscf, mrpt); no atomterm imports -- the
substrate objects arrive built.  No chkfiles, no caching: everything
derives from the arguments in one deterministic pass.
"""

import numpy as np
from pyscf import fci, mcscf, mrpt
from pyscf.mrpt import nevpt2 as _nevpt2

__all__ = [
    "ClassRow",
    "PT2Row",
    "readout_mc",
    "assert_stored_roots",
    "nevpt2_rows",
]

# the Koopmans classes in the kernel's own order, with its print labels
# and its per-class call signature (the kernel's argument shapes)
_CLASS_CALLS = (
    ("Sr", "Sr (-1)'", "cidms"),
    ("Si", "Si (+1)'", "cidms"),
    ("Sijrs", "Sijrs (0)", "eris"),
    ("Sijr", "Sijr (+1)", "dms"),
    ("Srsi", "Srsi (-1)", "dms"),
    ("Srs", "Srs (-2)", "dms"),
    ("Sij", "Sij (+2)", "dms"),
    ("Sir", "Sir (0)'", "dms"),
)

_VARIANTS = ("nevpt2-sc",)


class ClassRow:
    """One Koopmans class of the strongly contracted decomposition.

    ``label`` -- the kernel's own tag (``Sr``, ``Sijrs``, ...);
    ``note`` -- the kernel's print label (the charge suffix);
    ``norm`` -- the class's perturber norm (nonnegative);
    ``energy`` -- the class's E(2) contribution.
    """

    def __init__(self, label: str, note: str, norm: float, energy: float):
        self.label = label
        self.note = note
        self.norm = norm
        self.energy = energy

    @property
    def e_per_norm(self) -> float:
        """The compression meter: energy per unit perturber norm.

        SC's softened intruder read -- a class compressing a large
        coupling into a small energy (|ratio| far below its peers)
        sits near a small denominator; a zero norm carries no meter.
        """
        return self.energy / self.norm if self.norm else float("nan")

    def __repr__(self):
        return (f"ClassRow({self.label} {self.note}: "
                f"norm {self.norm:.6e}, e {self.energy:.6e})")


class PT2Row:
    """One root's PT2 result: the correction + the class breakdown.

    ``e_ref`` -- the root's CASCI total (the readout's own value);
    ``e_corr`` -- the class-summed E(2); ``e_total`` -- their sum.
    """

    def __init__(self, root: int, variant: str, e_ref: float,
                 e_corr: float, classes: list[ClassRow]):
        self.root = root
        self.variant = variant
        self.e_ref = e_ref
        self.e_corr = e_corr
        self.classes = classes

    @property
    def e_total(self) -> float:
        return self.e_ref + self.e_corr

    def __repr__(self):
        return (f"PT2Row(root {self.root}, {self.variant}: "
                f"e_ref {self.e_ref:.10f}, e_corr {self.e_corr:.10f})")


def _root_energies(mc) -> np.ndarray:
    """Total energies of all roots (active-space e_cas + the common
    core part) -- the store's own root-energy convention."""
    e_cas = np.atleast_1d(np.asarray(mc.e_cas, dtype=float)).ravel()
    e_tot = float(np.atleast_1d(np.asarray(mc.e_tot,
                                           dtype=float)).ravel()[0])
    return e_cas + (e_tot - float(e_cas[0]))


def readout_mc(mol, mf, mo_full, ncas: int, nelecas: tuple[int, int],
               nroots: int):
    """The multi-root CASCI readout on the supplied orbitals.

    ``mo_full``'s column layout is the partition (core prefix, active
    block at ``[ncore : ncore+ncas]``, virtual tail); the kernel runs
    exactly as the stored cells' readout pass did.

    Returns the finished ``mcscf.CASCI`` object (``.ci`` a list at
    ``nroots > 1``) -- the object ``mrpt.NEVPT`` consumes.
    """
    mc = mcscf.CASCI(mf, ncas, nelecas)
    mc.mo_coeff = mo_full
    mc.fcisolver.nroots = nroots
    mc.kernel()
    return mc


def _projector(cis) -> np.ndarray:
    """The root set's joint projector -- the gauge-invariant object."""
    dim = cis[0].shape[0]
    p = np.zeros((dim, dim))
    for c in cis:
        v = np.asarray(c).ravel()
        p += np.outer(v, v.conj())
    return p


def assert_stored_roots(mc, ref_energies=None, ref_ci=None,
                        e_tol: float = 1e-11, c_tol: float = 1e-8,
                        o_tol: float = 1e-3, p_tol: float = 1e-6) -> None:
    """Refuse unless the readout reproduced the stored solution.

    Energies compare strict (absolute tolerance; the fresh readout's
    observed deviation sits at 1e-15, far below it).  CI vectors
    compare per-root through three clauses, each a declared class:
    amplitude equality (c_tol); the degenerate block's lawful gauge
    freedom (the joint projector, its own band p_tol -- an SA draw's
    solver tails move the block projector at the tail scale, an
    order below the vectors' mixing angles); and the single-root
    re-solve's davidson noise (the re-solve ANGLE -- |sin| between
    the vectors, linear in the unconverged mixing -- within o_tol of
    zero; the state's meters and E(2) are flat across it, while any
    physical mixture fails the energy clause first).

    Parameters
    ----------
    e_tol, c_tol, o_tol, p_tol : the refusal bands (absolute); force
        them to pin the band behavior (the defaults are the
        production refusal).

    Raises
    ------
    ValueError : on any mismatch -- the caller's store and the
        readout disagree about what the solution is.
    """
    if ref_energies is not None:
        e = _root_energies(mc)
        ref = np.atleast_1d(np.asarray(ref_energies, dtype=float)).ravel()
        if e.shape != ref.shape or not np.allclose(e, ref, atol=e_tol,
                                                   rtol=0.0):
            shown = ("shape mismatch"
                     if e.shape != ref.shape
                     else f"max diff {np.max(np.abs(e - ref)):.3e}")
            raise ValueError(
                "stored-roots refusal: the readout energies "
                f"({shown}) do not "
                "reproduce the stored energies -- the substrate does "
                "not reconstruct the stored run; PT2 refused")
    if ref_ci is not None:
        new = [mc.ci] if mc.fcisolver.nroots == 1 else list(mc.ci)
        ref = ([ref_ci] if np.asarray(ref_ci).ndim == 1
               else [np.asarray(c) for c in ref_ci])
        if len(new) != len(ref):
            raise ValueError(
                f"stored-roots refusal: {len(new)} roots read out vs "
                f"{len(ref)} stored")
        # amplitudes are the assertion, not the container shape (the
        # store keeps flat vectors, the reader returns (n, 1) columns)
        def _amplitudes_match(c, r) -> bool:
            if c.size != r.size:
                return False
            c, r = c.ravel(), r.ravel()
            if np.allclose(c, r, atol=c_tol, rtol=0.0):
                return True
            # the phase/noise clause: same state up to a global phase
            # and sub-percent davidson noise -- the re-solve ANGLE
            # (|sin|, linear in the mixing) within o_tol of zero.
            # The sine comes from the orthogonal residual: at small
            # angles sqrt(1-cos^2) cancels catastrophically
            cn, rn = np.linalg.norm(c), np.linalg.norm(r)
            if not cn or not rn:
                return False
            u = c / cn
            rn_u = r / rn
            sino = float(np.linalg.norm(rn_u - np.vdot(u, rn_u) * u))
            return bool(sino <= o_tol)

        per_root = all(_amplitudes_match(c, r) for c, r in zip(new, ref))
        if per_root:
            return
        if (all(c.size == r.size for c, r in zip(new, ref))
                and np.allclose(_projector(new), _projector(ref),
                                atol=p_tol, rtol=0.0)):
            return
        raise ValueError(
            "stored-roots refusal: the readout CI vectors match "
            "neither the stored roots nor their joint projector -- "
            "the solution is not the stored one; PT2 refused")


def _root_row(mc, root: int, variant: str) -> PT2Row:
    """One root's row: the canonicalization prelude + the eight class
    solves (the kernel's own recipe, mirrored; the sum is pinned
    against ``NEVPT(...).kernel()`` in the tests)."""
    e_ref = float(_root_energies(mc)[root])
    pt = mrpt.NEVPT(mc, root=root, density_fit=False)
    if pt.fcisolver.nroots > 1:
        # detach the shared list: the kernel's own loop writes the
        # canonicalized vector into slot ``root``; the caller's ``mc``
        # must stay untouched for the other roots
        pt.ci = list(pt.ci)
    pt.mo_coeff, ci_vec, pt.mo_energy = pt.canonicalize(
        pt.mo_coeff, ci=pt.load_ci(), cas_natorb=True, verbose=pt.verbose)
    if pt.fcisolver.nroots == 1:
        pt.ci = ci_vec
    else:
        pt.ci[root] = ci_vec

    ncore, ncas = pt.ncore, pt.ncas
    nocc = ncore + ncas
    ci = pt.load_ci()
    dm1, dm2, dm3 = fci.rdm.make_dm123("FCI3pdm_kern_sf", ci, ci, ncas,
                                       pt.nelecas)
    dms = {"1": dm1, "2": dm2, "3": dm3, "4": None}
    link_a = fci.cistring.gen_linkstr_index(range(ncas), pt.nelecas[0])
    link_b = fci.cistring.gen_linkstr_index(range(ncas), pt.nelecas[1])
    eris = _nevpt2._ERIS(pt, pt.mo_coeff)
    aaaa = eris["ppaa"][ncore:nocc, ncore:nocc].copy()
    dms["f3ca"] = _nevpt2._contract4pdm("NEVPTkern_cedf_aedf", aaaa, ci,
                                        ncas, pt.nelecas, (link_a, link_b))
    dms["f3ac"] = _nevpt2._contract4pdm("NEVPTkern_aedf_ecdf", aaaa, ci,
                                        ncas, pt.nelecas, (link_a, link_b))
    classes = []
    for fname, label, kind in _CLASS_CALLS:
        fn = getattr(_nevpt2, fname)
        if kind == "cidms":
            norm, energy = fn(pt, ci, dms, eris)
        elif kind == "dms":
            norm, energy = fn(pt, dms, eris)
        else:
            norm, energy = fn(pt, eris)
        classes.append(ClassRow(fname, label, float(norm), float(energy)))
    return PT2Row(root=root, variant=variant, e_ref=e_ref,
                  e_corr=float(sum(c.energy for c in classes)),
                  classes=classes)


def nevpt2_rows(mol, mf, mo_full, ncas: int, nelecas: tuple[int, int],
                nroots: int, *, roots=None, ref_energies=None, ref_ci=None,
                variant: str = "nevpt2-sc") -> list[PT2Row]:
    """PT2 rows for the declared roots of one solution.

    Parameters
    ----------
    mol, mf, mo_full : the substrate (``mf`` a converged-pattern SCF
        object; ``mo_full`` the full column partition).
    ncas, nelecas, nroots : the CAS shape and the readout's root count
        (``nroots`` must reproduce the store's own readout when refs
        are given).
    roots : the declared root indices (default: all).
    ref_energies, ref_ci : the stored solution's energies / CI stack
        (PySCF order) -- given, the readout is asserted against them
        first (fail-fast); omitted, the solution is taken as fresh.
    variant : the row's variant tag (the whitelist grows with the
        constructors).

    Returns
    -------
    One :class:`PT2Row` per declared root, in the declared order.

    Raises
    ------
    ValueError : unknown variant, or the stored-roots refusal.
    """
    if variant not in _VARIANTS:
        raise ValueError(
            f"unknown PT2 variant {variant!r}; available: {_VARIANTS} "
            "(further families join as constructors at this seam)")
    mc = readout_mc(mol, mf, mo_full, ncas, nelecas, nroots)
    assert_stored_roots(mc, ref_energies=ref_energies, ref_ci=ref_ci)
    return [_root_row(mc, k, variant)
            for k in (range(nroots) if roots is None else roots)]
