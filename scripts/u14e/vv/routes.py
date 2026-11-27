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
the C1 reference: the independent anchor.  The PINNED routes are the
production wiring (marked: under test) -- one ``wfnsym``-pinned
solver per D2h component of the term, on the D2h-tagged wrapper.  The
free readout has the shape of PySCF's own recommendation for a
state-averaged wavefunction (``examples/mrpt/41-for_state_average.py``):
optimize, then a multi-root CASCI on the converged orbitals.

Layering
--------
:class:`State` / :class:`Solution` / :class:`PricedState` /
:class:`Priced` -- the vocabulary; :func:`assert_consistent` -- the
birth certificate; ``free_*`` / ``pinned_*`` -- the routes;
``_optimize`` / ``_readout`` / ``_born`` / ``_pinned_states`` -- the
two births.
"""

from dataclasses import dataclass, replace

import numpy as np
from pyscf import ao2mo, fci, mcscf
from pyscf.fci import cistring

from atomterm.casscf import irrep_solvers, run_irrep_casci, run_irrep_sa
from atomterm.core import closed_shell_core_energy
from atomterm.irreps import orb_irreps
from atomterm.operators import active_space_integrals

from vv.meters import Meters
from vv.pricing import price_classes, price_pyscf
from vv.table import Route, Row, conv_of, fmt_range

# the optimizer budget (the user's directive: "enable even 200 iters")
# -- 4x PySCF's default 50-macro cap; the deep draws need it (observed:
# 192 macro cycles at one (6,13) draw)
BUDGET = {"max_cycle_macro": 200}
ORTHO_TOL = 1e-8     # max|C^T S C - I| of the orbitals a solution carries
EIGENPAIR_TOL = 1e-7  # Eh; |<Psi|H|Psi> - e_cas| of a consistent state
                    # (observed on consistent states: ~1e-12)


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
    certificate's readings (:func:`assert_consistent`).
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
                return price_pyscf(st), None, *m.of(st.ci)
            row = pricer(st)
            return row.e_corr, tuple(row.classes), *m.of(st.ci)

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
    (None through PySCF's kernel), <<L^2>>, <<S^2>>."""

    state: State
    e2: float
    classes: tuple | None
    l2: float
    s2: float

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

    def row(self, cas: str) -> Row:
        """The table row, in the route's own convention: a free row
        carries the first priced root; a pinned row carries the mean
        over the term's component blocks (E_method: the optimizer's
        ensemble energy)."""
        sol = self.solution
        cells = (fmt_range(self.l2), f"{max(self.s2):.4f}")
        if sol.route is Route.FREE:
            first = self.states[0]
            return Row(sol.method, sol.route, cas, first.state.e_ref,
                       first.e2, *cells, sol.conv)
        e_method = (sol.ensemble if sol.ensemble is not None
                    else float(np.mean([s.state.e_ref for s in self.states])))
        return Row(sol.method, sol.route, cas, e_method,
                   float(np.mean(self.e2)), *cells, sol.conv)


_PRICERS = {Route.FREE: price_pyscf, Route.PINNED: price_classes,
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
    fold = "raw" if route is Route.FREE else "averaged"
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


def _readout(mf, space, mo_coeff, nroots, ci0=None):
    """PySCF's multi-root CASCI on the given orbitals -- the
    readout of an optimized state-averaged wavefunction
    (``examples/mrpt/41-for_state_average.py``), optionally seeded
    from the optimizer's own vectors.  ``canonicalization`` off: the
    window exactly as converged."""
    cas = mcscf.CASCI(mf, space.ncas, space.sector)
    cas.mo_coeff = mo_coeff
    cas.canonicalization = False
    cas.fcisolver.nroots = nroots
    cas.kernel(ci0=ci0)
    return cas


def free_casci(sub, space) -> Solution:
    """The free route at the CI level: PySCF's CASCI in the window as
    built -- the frozen pricing class (what the perturbation sees
    before any relaxation)."""
    space.certify(sub.mf)
    cas = mcscf.CASCI(sub.mf, space.ncas, space.sector)
    cas.canonicalization = False
    cas.kernel(space.mo)
    conv = conv_of(cas.fcisolver.conv_tol, None)   # a CI solve: the requested tolerance
    return _born(Route.FREE, "casci", space, [State(cas, 0)], conv)


def free_casscf(sub, space, budget=BUDGET) -> Solution:
    """The free route, single root: PySCF's plain CASSCF from the
    window, then the readout of its lowest root in the converged
    orbitals."""
    space.certify(sub.mf)
    mc, conv = _optimize(sub.mf, space, budget)
    cas = _readout(sub.mf, space, mc.mo_coeff, nroots=1)
    return _born(Route.FREE, "casscf", space, [State(cas, 0)], conv)


def free_sa(sub, space, nroots, budget=BUDGET) -> Solution:
    """The free route, state-averaged: PySCF's state average over the
    ``nroots`` lowest roots (uniform weights), then the multi-root
    readout in the converged orbitals, seeded from the optimizer's
    own vectors (unseeded, the davidson tails mix a quasi-degenerate
    block at the 1e-5 level)."""
    space.certify(sub.mf)
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
    ``readout_shift`` reading.
    """
    space.certify(sub.mf)
    mol, mf = sub.d2h
    out = run_irrep_sa(mol, mf, space.active, space.ncore,
                       *space.sector, *term, conv_tol=conv_tol,
                       c_core=space.core, mc_opts=budget)
    mo = np.asarray(out["mo"])
    active = slice(space.ncore, space.ncore + space.ncas)
    roots = run_irrep_casci(mol, mf, mo[:, active], space.ncore,
                            *space.sector, *term, conv_tol=conv_tol,
                            c_core=mo[:, :space.ncore])
    conv = conv_of(out["conv_dE"], out["cycles"], bool(out["converged"]))
    states = _pinned_states(sub, space, mo, term, roots, conv_tol)
    sol = _born(Route.PINNED, "sa-casscf", space, states, conv)
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
    core_part = closed_shell_core_energy(mf, space.core)
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
           "pinned_sa"]
