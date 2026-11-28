"""The routes: what runs on an active space and what comes out.

THE RULE.  A route solves, and returns states that carry the orbitals
their CI vectors were solved in.  A :class:`State` is assembled only
inside this module, from one solve; no caller hands in orbitals and
vectors separately -- so the inconsistent pair that once mispriced
six table rows (the CI vectors of the CONVERGED window placed on the
START window's columns, a non-orthonormal set) cannot be built.  On
birth every solution takes its certificate (:func:`assert_consistent`).

Two route families, deliberately different instruments.  The FREE
routes are PySCF-direct -- plain CASCI / CASSCF / state average on
the C1 reference, and ONLY there: on a symmetry-tagged substrate
PySCF's symmetry-adapted classes pin the FCI to one irrep, a
different calculation (the free routes refuse it).  They are the
independent anchor.  The PINNED routes are the production wiring
(marked: under test) -- one ``wfnsym``-pinned solver per D2h
component of the term, on the D2h-tagged wrapper.  The
free readout has the shape of PySCF's own recommendation for a
state-averaged wavefunction (``examples/mrpt/41-for_state_average.py``):
optimize, then a multi-root CASCI on the converged orbitals.

THE RADIAL ROUTE (:func:`radial_sa`) is the state average CONSTRAINED
to the window's composition: the production driver
(``atomterm.casscf.constrained_casscf``) on the window's lm-pure
sector basis, under the radial mask -- a rotation only within an
(l, member) sector, the core frozen.  The free state average does not
keep a composition: given empty non-f orbitals it rotates them into f
(observed: an s and a p channel of the U14e base converge with
f-weight 0.999..1.000, the term's roots splitting by 0.4..1.3 mEh).
The route is seeded and checked, because a state average of ``nroots``
LOWEST roots can track the wrong set:
  - a cold Davidson solve can MISS components of an exactly
    degenerate term -- its P-space guess does not cover every M_L
    block (observed in closed U14e windows: 14 of 17 with 17 roots,
    16 of 17 with 30); the cold solves here are outright
    diagonalizations up to ``FULL_DIM`` determinants, and the seed is
    read wider than the term (:func:`term_seed`);
  - a warm-started solve keeps tracking a set it is given even when
    the window holds lower states (observed: a 5f^2 6d^2 term tracked
    30 mEh above lower states of its own window); at convergence a
    fresh wider solve must show the tracked roots as the window's
    lowest level, one degenerate term.

Layering
--------
:class:`State` / :class:`Solution` / :class:`PricedState` /
:class:`Priced` -- the vocabulary; :func:`assert_consistent` -- the
birth certificate; ``free_*`` / ``pinned_*`` / :func:`radial_sa` --
the routes; ``_optimize`` / ``_readout`` / ``_born`` /
``_pinned_states`` / ``_cold`` -- the births and the seeds'
instrument; ``_lowest_level`` -- the trap check.
"""

from dataclasses import dataclass, replace

import numpy as np
from pyscf import ao2mo, fci, mcscf
from pyscf.fci import cistring

from atomterm.casscf import (constrained_casscf, irrep_solvers,
                             run_irrep_casci, run_irrep_sa)
from atomterm.core import closed_shell_core_energy
from atomterm.detbasis import sector_dimension
from atomterm.irreps import orb_irreps
from atomterm.operators import active_space_integrals, window_closure
from atomterm.windows import mask_radial_rotations

from vv.meters import Meters
from vv.pricing import price_classes, price_pyscf, resorted
from vv.windows import ActiveSpace, lm_basis, orthonormalize
from vv.table import Route, Row, conv_of, fmt_range

# the optimizer budget (the user's directive: "enable even 200 iters")
# -- 4x PySCF's default 50-macro cap; the deep draws need it (observed:
# 192 macro cycles at one (6,13) draw)
BUDGET = {"max_cycle_macro": 200}
ORTHO_TOL = 1e-8     # max|C^T S C - I| of the orbitals a solution carries
EIGENPAIR_TOL = 1e-7  # Eh; |<Psi|H|Psi> - e_cas| of a consistent state
                    # (observed on consistent states: ~1e-12)
SEED_MARGIN = 13     # the extra roots of the cold solve a state average
                     # is seeded from -- enough to see the window's next
                     # levels beside the term
SEED_DEV = 1e-4      # a seed root's distance to the term's (<L^2>, <S^2>)
LEVEL_TOL = 1e-3     # mEh: the spread of the tracked term's roots in the
                     # fresh solve (a closed window's term: <= 2e-4)
TRAP_TOL = 1e-6      # Eh: the tracked ensemble against the fresh solve's
                     # lowest level (observed <= 2e-9 when they agree)
FULL_DIM = 20000    # determinants: up to this sector size a cold solve is
                     # PySCF's outright diagonalization (every determinant
                     # in its P-space; ncas 27 in the (4, 0) sector) --
                     # above it a Davidson again, which may miss a
                     # component of an exactly degenerate term


# ------------------------------------------------------------ the vocabulary

@dataclass(frozen=True)
class State:
    """One CI root on the orbitals it was solved in -- the unit the
    pricing and the meters consume.

    ``mc`` is NEVPT's own input: a kernel-ready CASCI object carrying
    the full MO matrix (the partition layout [core | active |
    virtual]), the CI roots and the CAS energies.  ``root`` selects
    one of them; ``irrep`` names the pinned block's D2h component
    (None on the free routes).
    """

    mc: object
    root: int = 0
    irrep: str | None = None

    @property
    def mol(self):
        return self.mc._scf.mol

    @property
    def mf(self):
        return self.mc._scf

    @property
    def mo_coeff(self) -> np.ndarray:
        return np.asarray(self.mc.mo_coeff)

    @property
    def ncore(self) -> int:
        return int(self.mc.ncore)

    @property
    def ncas(self) -> int:
        return int(self.mc.ncas)

    @property
    def nelecas(self) -> tuple[int, int]:
        return tuple(int(x) for x in self.mc.nelecas)

    @property
    def active(self) -> np.ndarray:
        """The active block -- the columns this root was solved in."""
        return self.mo_coeff[:, self.ncore:self.ncore + self.ncas]

    @property
    def core(self) -> np.ndarray:
        return self.mo_coeff[:, :self.ncore]

    @property
    def ci(self) -> np.ndarray:
        """The root's CI vector, PySCF order and shape (nsa, nsb)."""
        c = self.mc.ci
        if isinstance(c, np.ndarray) and c.ndim == 2:   # nroots == 1
            return c
        return np.asarray(c[self.root])

    @property
    def e_cas(self) -> float:
        return float(np.atleast_1d(np.asarray(self.mc.e_cas)).ravel()[self.root])

    @property
    def e_ref(self) -> float:
        """The root's total energy: ``e_cas`` plus the solution's core
        bookkeeping (the first root's total minus its CAS energy)."""
        cas = np.atleast_1d(np.asarray(self.mc.e_cas)).ravel()
        tot0 = float(np.atleast_1d(np.asarray(self.mc.e_tot)).ravel()[0])
        return self.e_cas + (tot0 - float(cas[0]))


@dataclass(frozen=True)
class Solution:
    """What one route returns.

    ``states`` -- the roots in the readout's order (on the free routes
    the single, or the lowest, root is the table's); ``conv`` -- the
    conv cell; ``ensemble`` -- the optimizer's ensemble energy on the
    state-averaged routes, else None; ``certificate`` -- the birth
    certificate's readings (:func:`assert_consistent`, plus the
    route's own: ``core_leak``; on the pinned SA also ``ortho_in``,
    ``readout_shift``, ``ensemble_dev``).
    """

    route: Route
    method: str
    states: tuple[State, ...]
    conv: str
    ensemble: float | None = None
    certificate: dict | None = None

    @property
    def nroots(self) -> int:
        return len(self.states)

    def meters(self, roots=None) -> tuple[tuple[float, float], ...]:
        """(<<L^2>>, <<S^2>>) per declared root (default: all), one
        meter per distinct orbital set."""
        return tuple(self._per_state(roots, lambda st, m: m.of(st.ci)))

    def price(self, roots=None, how: str | None = None) -> "Priced":
        """Price and meter the declared roots (default: all).

        Parameters
        ----------
        roots : the declared root indices (default: all).
        how : ``"pyscf"`` -- PySCF's own kernel (the free routes'
            pricer, the external anchor); ``"classes"`` -- the
            engine's recipe with the Koopmans class breakdown (the
            pinned route's pricer).  Default: the route's own.
        """
        pricer = _PRICERS[how] if how else _PRICERS[self.route]

        def one(st, m):
            if pricer is price_pyscf:
                return price_pyscf(st), None, *m.of(st.ci), resorted(st)
            row = pricer(st)
            return row.e_corr, tuple(row.classes), *m.of(st.ci), resorted(st)

        return Priced(self, tuple(PricedState(st, *val)
                                  for st, val in zip(self._selected(roots),
                                                     self._per_state(roots, one))))

    def _selected(self, roots):
        idx = range(self.nroots) if roots is None else roots
        return [self.states[k] for k in idx]

    def _per_state(self, roots, f):
        meters, out = {}, []
        for st in self._selected(roots):
            m = meters.get(id(st.mo_coeff))
            if m is None:
                m = meters[id(st.mo_coeff)] = Meters(st.mol, st.active, st.nelecas)
            out.append(f(st, m))
        return out


@dataclass(frozen=True)
class PricedState:
    """One priced, metered root: E(2) in Eh, the Koopmans class rows
    (None through PySCF's kernel), <<L^2>>, <<S^2>>, and ``resorted``
    -- how many empty active orbitals the pricing prelude re-partitions
    with the virtuals (:func:`vv.pricing.resorted`; 0 when the priced
    partition is the solved window)."""

    state: State
    e2: float
    classes: tuple | None
    l2: float
    s2: float
    resorted: int = 0

    @property
    def e_total(self) -> float:
        return self.state.e_ref + self.e2


@dataclass(frozen=True)
class Priced:
    """A solution with its declared roots priced and metered."""

    solution: Solution
    states: tuple[PricedState, ...]

    @property
    def e2(self) -> list[float]:
        return [s.e2 for s in self.states]

    @property
    def l2(self) -> list[float]:
        return [s.l2 for s in self.states]

    @property
    def s2(self) -> list[float]:
        return [s.s2 for s in self.states]

    @property
    def resorted(self) -> int:
        """The largest :attr:`PricedState.resorted` over the states (0:
        every priced partition is the solved window)."""
        return max((s.resorted for s in self.states), default=0)

    def row(self, cas: str) -> Row:
        """The table row, in the route's own convention: a free row
        carries the first priced root; a pinned row carries the mean
        over the term's component blocks (E_method: the optimizer's
        ensemble energy)."""
        sol = self.solution
        cells = (fmt_range(self.l2), f"{max(self.s2):.4f}")
        if sol.route is not Route.PINNED:
            first = self.states[0]
            return Row(sol.method, sol.route, cas, first.state.e_ref,
                       first.e2, *cells, sol.conv)
        e_method = (sol.ensemble if sol.ensemble is not None
                    else float(np.mean([s.state.e_ref for s in self.states])))
        return Row(sol.method, sol.route, cas, e_method,
                   float(np.mean(self.e2)), *cells, sol.conv)


_PRICERS = {Route.FREE: price_pyscf, Route.PINNED: price_classes,
            Route.RADIAL: price_pyscf,
            "pyscf": price_pyscf, "classes": price_classes}


# -------------------------------------------------------- the birth certificate

def core_leak(space, states) -> float:
    """How far the solution's core block drifted off the active
    space's frozen core span: the S-norm of the component outside the
    span (0 for a space used as built).  The optimizer returns its
    core block off the frozen span at its convergence noise (observed
    1.1e-06 on carbon, with and without canonicalization)."""
    worst = 0.0
    for st in states:
        s = st.mf.get_ovlp()
        core = st.core
        d = core - space.core @ (space.core.T @ s @ core)
        worst = max(worst, float(np.sqrt(np.trace(d.T @ s @ d))))
    return worst

FOLDS = ("raw", "averaged")


def cas_hamiltonian(mol, mf, active, core, fold: str = "raw") -> tuple:
    """The active-space Hamiltonian: (h1, eri) in the given active
    columns with the core folded into h1 -- PySCF's own FCI feed.

    Parameters
    ----------
    fold : ``"raw"`` -- the core field as the SCF solved it (PySCF's
        CASSCF/CASCI convention; the free routes' solvers);
        ``"averaged"`` -- the core density's l = 0 component
        (``operators.active_space_integrals``'s default; the pinned
        solvers').  The two are the SAME Hamiltonian on a purified
        core (spherical by construction) and differ by the core's
        polarization energy on a raw one (observed on a frozen
        carbon ROHF core: 2.9e-05 Eh).
    """
    h1, _ = active_space_integrals(mol, mf, active, core.shape[1],
                                   spherical_core=(fold == "averaged"),
                                   c_core=core)
    return h1, ao2mo.kernel(mol, active)


def energy_expectation(mol, mf, active, core, sector, ci,
                       fold: str = "raw") -> float:
    """<<Psi|H|Psi>> in the given active block -- the CAS energy's
    convention.  For a state's own columns this is its ``e_cas``;
    on any other columns it is strictly higher (the variational
    principle is the detector of a foreign vector)."""
    h1, eri = cas_hamiltonian(mol, mf, active, core, fold)
    return float(fci.direct_spin1.energy(h1, eri, ci, active.shape[1], sector))


def assert_consistent(states, space=None, fold: str = "raw",
                      tol=EIGENPAIR_TOL) -> dict:
    """The birth certificate: refuse unless every state is a
    CONSISTENT object -- what NEVPT actually consumes.

    (1) its orbitals are S-orthonormal (``ORTHO_TOL``) and the CASCI
        object's partition (ncore, ncas, nelecas) is the active
        space's, when one is given;
    (2) its CI vector is stationary in ITS OWN orbitals: the energy
        expectation in the state's active block equals the stored CAS
        energy within ``tol``.

    Parameters
    ----------
    states : the solution's states.
    space : the :class:`~vv.windows.ActiveSpace` the routes ran on.
    fold : the core fold of the SOLVERS that produced the states
        (:func:`cas_hamiltonian`) -- each route passes its own.
    tol : the energy residual's band, Eh.

    Returns
    -------
    ``{"ortho": the worst orthonormality deviation, "residual": the
    worst energy residual}``.

    Raises
    ------
    ValueError : a state fails either clause; the message names the
        state and the reading.
    """
    ortho, residual, hamiltonians = 0.0, 0.0, {}
    for st in states:
        mo = st.mo_coeff
        if space is not None and (st.ncore, st.ncas, st.nelecas) != (
                space.ncore, space.ncas, tuple(space.sector)):
            raise ValueError(
                f"state {st.root}: its CASCI object's partition "
                f"({st.ncore}, {st.ncas}, {st.nelecas}) is not the "
                f"active space's ({space.ncore}, {space.ncas}, "
                f"{tuple(space.sector)})")
        s = st.mf.get_ovlp()
        dev = float(np.abs(mo.T @ s @ mo - np.eye(mo.shape[1])).max())
        ortho = max(ortho, dev)
        if dev > ORTHO_TOL:
            raise ValueError(
                f"state {st.root}: its orbitals are not S-orthonormal "
                f"(max|C^T S C - I| = {dev:.2e} > {ORTHO_TOL:g})")
        key = id(st.mc.mo_coeff)
        if key not in hamiltonians:
            hamiltonians[key] = cas_hamiltonian(st.mol, st.mf, st.active,
                                                st.core, fold)
        h1, eri = hamiltonians[key]
        e = float(fci.direct_spin1.energy(h1, eri, st.ci, st.ncas, st.nelecas))
        res = abs(e - st.e_cas)
        residual = max(residual, res)
        if res > tol:
            raise ValueError(
                f"state {st.root}: its CI vector is not stationary in "
                f"its own orbitals -- <Psi|H|Psi> misses the stored "
                f"CAS energy by {res:.2e} Eh (band {tol:g}): the "
                f"vector belongs to another orbital set")
    return {"ortho": ortho, "residual": residual}


# ----------------------------------------------------------------- the births

def _born(route, method, space, states, conv, ensemble=None) -> Solution:
    # each route's states are certified in its own solvers' core fold
    fold = "averaged" if route is Route.PINNED else "raw"
    cert = assert_consistent(states, space, fold)
    cert["core_leak"] = core_leak(space, states)
    return Solution(route, method, tuple(states), conv, ensemble, cert)


def _optimize(mf, space, budget, weights=None) -> tuple:
    """PySCF's CASSCF -- state-averaged over ``weights`` when given --
    run to its budget (``conv_tol`` stays PySCF's 1e-7).  Returns
    (mc, the conv cell): the cell's |dE| is the LAST callback's (the
    final macro-cycle end) and the count is the TRUE macro-cycle
    count (PySCF's ``imacro`` is 1-based and its callback fires per
    MICRO iteration AND per macro cycle)."""
    mc = mcscf.CASSCF(mf, space.ncas, space.sector)
    for k, val in budget.items():
        setattr(mc, k, val)
    mc.conv_tol = 1e-7
    trace = {"imacro": 0, "de": float("nan")}

    def _cb(envs):
        trace["de"] = abs(float(envs.get("de", float("nan"))))
        m = int(envs.get("imacro", 0))
        if m > trace["imacro"]:
            trace["imacro"] = m

    mc.callback = _cb
    if weights is not None:
        mc.state_average_(list(weights))
    mc.kernel(space.mo)
    return mc, conv_of(trace["de"], trace["imacro"], bool(mc.converged))


def _readout(mf, space, mo_coeff, nroots, ci0=None, exact=False):
    """PySCF's multi-root CASCI on the given orbitals -- the
    readout of an optimized state-averaged wavefunction
    (``examples/mrpt/41-for_state_average.py``), optionally seeded
    from the optimizer's own vectors.  ``canonicalization`` off: the
    window exactly as converged.

    ``exact`` -- an outright diagonalization instead of the Davidson:
    PySCF diagonalizes its P-space Hamiltonian when the P-space holds
    every determinant of the sector, so the solver's ``pspace_size``
    is set to the sector's dimension (no guess, no missed root)."""
    cas = mcscf.CASCI(mf, space.ncas, space.sector)
    cas.mo_coeff = mo_coeff
    cas.canonicalization = False
    cas.fcisolver.nroots = nroots
    if exact:
        cas.fcisolver.pspace_size = sector_dimension(space.ncas, *space.sector)
    cas.kernel(ci0=ci0)
    return cas


def _certify_free(sub, space) -> None:
    """The free and the radial routes' substrate is the C1 reference.
    On a symmetry-tagged molecule PySCF hands a plain CASCI/CASSCF its
    symmetry-adapted classes, and the FCI is pinned to one irrep --
    a different calculation (observed 2e-02 Eh above the sector
    ground on one window); the symmetry-tagged substrate belongs to
    the pinned routes.

    Raises
    ------
    ValueError : the substrate's molecule is symmetry-tagged; or the
        active space fails its certificate.
    """
    if getattr(sub.mol, "symmetry", False):
        raise ValueError(
            "a free route on a symmetry-tagged substrate: PySCF pins "
            "the FCI to one irrep there (a different calculation) -- "
            "the free routes run on the C1 reference")
    space.certify(sub.mf)


def free_casci(sub, space) -> Solution:
    """The free route at the CI level: PySCF's CASCI in the window as
    built -- the frozen pricing class (what the perturbation sees
    before any relaxation)."""
    _certify_free(sub, space)
    cas = mcscf.CASCI(sub.mf, space.ncas, space.sector)
    cas.canonicalization = False
    cas.kernel(space.mo)
    conv = conv_of(cas.fcisolver.conv_tol, None)   # a CI solve: the requested tolerance
    return _born(Route.FREE, "casci", space, [State(cas, 0)], conv)


def free_casscf(sub, space, budget=BUDGET) -> Solution:
    """The free route, single root: PySCF's plain CASSCF from the
    window, then the readout of its lowest root in the converged
    orbitals."""
    _certify_free(sub, space)
    mc, conv = _optimize(sub.mf, space, budget)
    cas = _readout(sub.mf, space, mc.mo_coeff, nroots=1)
    return _born(Route.FREE, "casscf", space, [State(cas, 0)], conv)


def _cold(sub, space, mo, nroots) -> tuple[list, list]:
    """A cold multi-root readout of ``min(nroots, the sector's
    dimension)`` roots in the window whose [core | active] is
    ``mo``'s: (the roots as uncategorized states, their (<L^2>,
    <S^2>)) -- the seeds' and the trap check's instrument.  Exact up
    to ``FULL_DIM`` determinants (the module docstring's first
    warning)."""
    dim = sector_dimension(space.ncas, *space.sector)
    nroots = min(nroots, dim)
    mo = orthonormalize(mo, sub.mf.get_ovlp(), space.ncore, space.ncas)
    win = ActiveSpace(mo, space.ncore, space.ncas, space.sector)
    cas = _readout(sub.mf, win, mo, nroots, exact=dim <= FULL_DIM)
    states = [State(cas, k) for k in range(nroots)]
    m = Meters(sub.mol, win.active, win.sector)
    return states, [m.of(st.ci) for st in states]


def term_seed(sub, space, term, nroots, margin=SEED_MARGIN) -> list:
    """The term's components in the window as built: from a cold
    solve of ``nroots + margin`` roots, the ``nroots`` lowest that
    carry the term (within ``SEED_DEV``) -- the starting vectors of a
    state average that must describe the term.

    Raises
    ------
    ValueError : fewer than ``nroots`` of the solved roots carry the
        term.
    """
    states, meters = _cold(sub, space, space.mo, nroots + margin)
    hit = [k for k, (l2, s2) in enumerate(meters)
           if max(abs(l2 - term[0]), abs(s2 - term[1])) <= SEED_DEV]
    if len(hit) < nroots:
        raise ValueError(
            f"the frozen window holds {len(hit)} of the term's "
            f"{nroots} components among its lowest {len(states)} "
            f"roots (<L^2> {[round(x[0], 4) for x in meters]})")
    return [states[k].ci for k in hit[:nroots]]


def _lowest_level(sub, space, mo, term, nroots, tracked) -> float:
    """The TRAP CHECK of a state average: in a fresh, wider cold solve
    of the window whose [core | active] is ``mo``'s, the lowest
    ``nroots`` roots must be ONE degenerate term (``LEVEL_TOL``,
    ``SEED_DEV``) lying at the tracked ensemble energy (``TRAP_TOL``).

    Returns
    -------
    The gap (mEh) between that level and the window's next one; NaN
    when the sector holds no level above it.

    Raises
    ------
    ValueError : any clause fails.
    """
    states, meters = _cold(sub, space, mo, nroots + SEED_MARGIN)
    e = [st.e_ref for st in states[:nroots]]
    spread = (max(e) - min(e)) * 1e3
    dev = max(max(abs(l2 - term[0]), abs(s2 - term[1]))
              for l2, s2 in meters[:nroots])
    if (spread > LEVEL_TOL or dev > SEED_DEV
            or abs(float(np.mean(e)) - tracked) > TRAP_TOL):
        raise ValueError(
            "the state average's tracked roots are not the converged "
            f"window's lowest {nroots}: a fresh solve reads a lowest "
            f"level of spread {spread:.4f} mEh (band {LEVEL_TOL:g}) at "
            f"max term deviation {dev:.1e} (band {SEED_DEV:g}), "
            f"{float(np.mean(e)) - tracked:+.2e} Eh from the tracked "
            f"ensemble (band {TRAP_TOL:g})")
    if len(states) <= nroots:
        return float("nan")
    return (states[nroots].e_ref - max(e)) * 1e3


def radial_sa(sub, space, term, nroots, budget=BUDGET) -> Solution:
    """The state average CONSTRAINED to the window's composition (the
    module docstring's radial route): the production driver on the
    lm-pure sector basis under the radial mask, seeded with the term's
    components, its readout seeded from the optimizer's own vectors.

    Parameters
    ----------
    term : (L(L+1), S(S+1)) -- the expectations every root carries.

    Returns
    -------
    The solution (``Route.RADIAL``); its certificate carries the
    birth readings, ``level_gap`` (:func:`_lowest_level`) and
    ``closure`` (the converged window's largest closure meter).

    Raises
    ------
    ValueError : the window is not lm-pure; the seed is refused; or
        the trap check (:func:`_lowest_level`) fails -- in a fresh
        solve of the converged window the tracked roots are not its
        lowest ``nroots``, one degenerate term.
    """
    _certify_free(sub, space)
    mol, mf = sub.mol, sub.mf
    mo, labels = lm_basis(sub, space)
    mask = mask_radial_rotations(labels, space.ncore, space.ncas)
    seed = term_seed(sub, space, term, nroots)
    out = constrained_casscf(mf, space.ncas, space.sector, None, mo_full=mo,
                             mask=mask, weights=[1.0 / nroots] * nroots,
                             mc_opts={**budget, "ci": seed})
    cas = _readout(mf, space, np.asarray(out.mo_coeff), nroots,
                   ci0=list(out.ci))
    conv = conv_of(out.conv_dE, out.cycles, bool(out.mc.converged))
    sol = replace(_born(Route.RADIAL, "sa-casscf", space,
                        [State(cas, k) for k in range(nroots)], conv),
                  ensemble=float(out.mc.e_tot))
    sol.certificate["level_gap"] = _lowest_level(
        sub, space, np.asarray(out.mo_coeff), term, nroots,
        float(np.mean([st.e_ref for st in sol.states])))
    sol.certificate["closure"] = max(window_closure(
        mol, sol.states[0].active).values())
    return sol


def free_sa(sub, space, nroots, budget=BUDGET) -> Solution:
    """The free route, state-averaged: PySCF's state average over the
    ``nroots`` lowest roots (uniform weights), then the multi-root
    readout in the converged orbitals, seeded from the optimizer's
    own vectors (unseeded, the davidson tails mix a quasi-degenerate
    block at the 1e-5 level)."""
    _certify_free(sub, space)
    mc, conv = _optimize(sub.mf, space, budget,
                         weights=[1.0 / nroots] * nroots)
    nsa = cistring.num_strings(space.ncas, space.sector[0])
    nsb = cistring.num_strings(space.ncas, space.sector[1])
    ci0 = [np.asarray(c).ravel().reshape(nsa, nsb) for c in mc.ci]
    cas = _readout(sub.mf, space, mc.mo_coeff, nroots, ci0)
    sol = _born(Route.FREE, "sa-casscf", space,
                [State(cas, k) for k in range(nroots)], conv)
    return replace(sol, ensemble=float(mc.e_tot))


def pinned_casci(sub, space, term, conv_tol: float = 1e-12) -> Solution:
    """The pinned route at the CI level: one block solve per D2h
    component of the term, in the window as built (the production
    runner ``run_irrep_casci`` -- marked: under test).

    Parameters
    ----------
    term : (L, parity) the target term (``parity`` is the
        configuration's total inversion parity).
    conv_tol : the pinned solvers' CI tolerance.
    """
    space.certify(sub.mf)
    roots = run_irrep_casci(*sub.d2h, space.active, space.ncore,
                            *space.sector, *term, conv_tol=conv_tol,
                            c_core=space.core)
    states = _pinned_states(sub, space, space.mo, term, roots, conv_tol)
    return _born(Route.PINNED, "casci", space, states,
                 conv_of(conv_tol, None))


def pinned_sa(sub, space, term, budget=BUDGET, conv_tol: float = 1e-12) -> Solution:
    """The pinned route, state-averaged: the component ensemble
    optimized (the production runner ``run_irrep_sa`` -- marked: under
    test), then the exact per-component roots re-solved in the
    CONVERGED window ON THE ORBITALS AS RETURNED -- the readout of
    PySCF's own recommendation for a state-averaged wavefunction
    (``examples/mrpt/41-for_state_average.py``).  The states are born
    from that readout: the optimizer's orbitals with their own core
    block, never the start window's.

    Notes
    -----
    The optimizer returns its core block off the frozen span at its
    convergence noise, and ``run_irrep_sa``'s internal readout folds
    the FROZEN core while using the returned active block -- a
    consistent-to-1e-6 object (:func:`core_leak`).  The states here
    are consistent exactly: their readout folds the returned core
    block.  The shift between the two readouts is the certificate's
    ``readout_shift`` reading; ``ortho_in`` is the deviation of the
    optimizer's orbitals as returned; ``ensemble_dev`` is the distance
    from the readout's component mean to the optimizer's ensemble
    energy (the two readouts describe one solution).
    """
    space.certify(sub.mf)
    mol, mf = sub.d2h
    out = run_irrep_sa(mol, mf, space.active, space.ncore,
                       *space.sector, *term, conv_tol=conv_tol,
                       c_core=space.core, mc_opts=budget)
    # the optimizer's orbitals, made exactly S-orthonormal block by
    # block (its symmetry-adapted optimizer returns the virtual block
    # orthonormal only to ~4e-05; the C1 one: 5e-12)
    returned = np.asarray(out["mo"])
    ortho_in = float(np.abs(returned.T @ mf.get_ovlp() @ returned
                            - np.eye(returned.shape[1])).max())
    mo = orthonormalize(returned, mf.get_ovlp(), space.ncore, space.ncas)
    active = slice(space.ncore, space.ncore + space.ncas)
    roots = run_irrep_casci(mol, mf, mo[:, active], space.ncore,
                            *space.sector, *term, conv_tol=conv_tol,
                            c_core=mo[:, :space.ncore])
    conv = conv_of(out["conv_dE"], out["cycles"], bool(out["converged"]))
    states = _pinned_states(sub, space, mo, term, roots, conv_tol)
    sol = _born(Route.PINNED, "sa-casscf", space, states, conv)
    sol.certificate["ortho_in"] = ortho_in
    # the readout's component mean against the optimizer's own ensemble
    sol.certificate["ensemble_dev"] = abs(float(np.mean(
        [s.e_ref for s in states])) - float(out["e_tot"]))
    sol.certificate["readout_shift"] = max(
        float(np.max(np.abs(np.asarray(e) - np.asarray(out["roots"][irrep][0]))))
        for irrep, (e, _ci) in roots.items())
    return replace(sol, ensemble=float(out["e_tot"]))


def _pinned_states(sub, space, mo, term, roots, conv_tol) -> list:
    """The states of a pinned readout: per component block a
    kernel-ready CASCI object on the ONE MO matrix the block's roots
    were solved in (``mo``), carrying the block's own solver (mrpt's
    canonicalization then stays in the block) and the closed-core
    bookkeeping in ``e_tot``."""
    mol, mf = sub.d2h
    mo = np.asarray(mo, dtype=float)
    orbsyms = orb_irreps(mol, mo[:, space.ncore:space.ncore + space.ncas])
    # the core bookkeeping of the core the readout folded: the states'
    # own core block
    core_part = closed_shell_core_energy(mf, mo[:, :space.ncore])
    states = []
    for irrep, solver, mult in irrep_solvers(mol, orbsyms, *term, conv_tol):
        e_arr, ci = roots[irrep]
        e_arr = np.atleast_1d(np.asarray(e_arr, dtype=float))
        mc = mcscf.CASCI(mf, space.ncas, space.sector)
        mc.mo_coeff = mo
        mc.fcisolver = solver
        mc.orbsym = orbsyms
        mc.fcisolver.orbsym = np.asarray(orbsyms)
        # PySCF's native form: a single-root object holds the bare
        # (nsa, nsb) vector (mrpt's load_ci returns mc.ci as it is
        # when nroots == 1), a multi-root one a list of them
        ci_arr = np.asarray(ci)
        mc.ci = (ci_arr[0] if ci_arr.ndim == 3 else ci_arr) if mult == 1 \
            else [np.asarray(x) for x in ci_arr]
        mc.e_cas = e_arr[:mult]
        mc.e_tot = float(mc.e_cas[0] + core_part)
        states += [State(mc, k, irrep) for k in range(mult)]
    return states


__all__ = ["BUDGET", "Priced", "PricedState", "Solution", "State",
           "assert_consistent", "cas_hamiltonian", "energy_expectation",
           "free_casci", "free_casscf", "free_sa", "pinned_casci",
           "pinned_sa", "radial_sa", "term_seed"]
