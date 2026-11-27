"""The routes: the free routes' two kernels (PySCF-direct), the pinned
routes (the production wiring, marked), the pinned solvers'
kernel-ready block objects, the lane route's substrate."""

import numpy as np  # noqa: E402
from pyscf import mcscf, mrpt  # noqa: E402

from vv.meters import meters  # noqa: E402

BUDGET_200 = {"max_cycle_macro": 200}  # the optimizer budget: the
# user's directive ("enable even 200 iters") — 4x PySCF's default
# 50-macro cap; the campaign's effective cap (the stored callback
# counts' '200' = 50 macro + ~150 micro) never converged the deep
# draws ((6,13): 192 macro needed)


# the derivations' thin helpers (one home; the legs show the setup) ----

def free_kernel(mf, ncas: int, sector, mo_full, weights=None,
                nroots: int = 1, track_iters: bool = False,
                max_cycle_macro: int = 200):
    """The free route's two kernels, written once: CASSCF (plain
    single root, or ``state_average_`` over ``weights``), then the
    two-pass readout (CASCI, ``canonicalization = False`` -- the
    window exactly as converged), seeded from the SA kernel's own
    vectors (the faithful two-pass form; unseeded, the davidson tails
    mix a quasi-degenerate block at the 1e-5 level and the meters
    carry it).  Returns (mc, cas).

    The optimizer budget is the user's directive ("enable even 200
    iters"): ``max_cycle_macro=200`` -- above the campaign cells'
    effective cap (PySCF's default 50 macro cycles; the campaign set
    no optimizer knobs, its stored ``cycles`` counts callback
    invocations = macro + micro).  ``conv_tol`` stays PySCF's 1e-7 --
    the campaign's effective value (converged draws landed |dE|
    1e-8..1e-9, below the tolerance, on the final macro step).

    ``track_iters`` rides the optimizer's callback and sets
    ``mc.n_macro``/``mc.final_de`` -- the conv record.  PySCF's
    ``imacro`` is 1-based and its callback fires per MICRO iteration
    AND per macro cycle, so: ``n_macro`` = the max ``imacro`` seen
    (the true macro-cycle count), ``final_de`` = the LAST call's
    |dE| (the last callback is always the final macro-cycle-end
    call; taking the first call of the final cycle instead reads the
    previous cycle's stale delta, and ``max imacro + 1`` overcounts).
    The default off leaves the historical path untouched."""
    mc = mcscf.CASSCF(mf, ncas, sector)
    mc.max_cycle_macro = max_cycle_macro
    mc.conv_tol = 1e-7
    if track_iters:
        _it = {"imacro": 0, "de": float("nan")}

        def _cb(envs):
            _it["de"] = abs(float(envs.get("de", float("nan"))))
            m = int(envs.get("imacro", 0))
            if m > _it["imacro"]:
                _it["imacro"] = m
        mc.callback = _cb
    if weights is not None and len(weights) > 1:
        mc.state_average_(list(weights))
    mc.kernel(mo_full)
    if track_iters:
        mc.n_macro = _it["imacro"]
        mc.final_de = _it["de"]
    cas = mcscf.CASCI(mf, ncas, sector)
    cas.mo_coeff = mc.mo_coeff
    cas.canonicalization = False
    cas.fcisolver.nroots = nroots
    if weights is not None and len(weights) > 1:
        from pyscf.fci import cistring
        nsa = cistring.num_strings(ncas, sector[0])
        nsb = cistring.num_strings(ncas, sector[1])
        cas.kernel(ci0=[np.asarray(c).ravel().reshape(nsa, nsb)
                        for c in mc.ci])
    else:
        cas.kernel()
    return mc, cas


def core_part_of(cas) -> tuple:
    """The per-root totals convention of the analyzation seam:
    ``e_cas[k] + core part`` (e_tot's first entry minus e_cas's
    first).  Returns (core_e, e_cas array)."""
    e_tot_0 = float(np.atleast_1d(np.asarray(cas.e_tot)).ravel()[0])
    e_cas = np.atleast_1d(np.asarray(cas.e_cas)).ravel()
    return e_tot_0 - float(e_cas[0]), e_cas


def price(cas, expect, roots=None) -> list:
    """PT2 through the engine's own NEVPT per root + the meters on
    the priced states (the readout's own roots)."""
    out = []
    e_cas = np.atleast_1d(np.asarray(cas.e_cas)).ravel()
    multi = len(e_cas) > 1
    for k in (range(len(e_cas)) if roots is None else roots):
        pt = mrpt.NEVPT(cas, root=k)
        if multi:
            pt.ci = list(pt.ci)   # detach: the kernel writes the shared list
        e2 = float(pt.kernel())
        # nroots=1 returns the (nsa, nsb) vector itself, not a list
        ci = np.asarray(cas.ci if not multi else cas.ci[k]).ravel()
        out.append({"root": k, "e2": e2, "l2": expect("l2", ci),
                    "s2": expect("s2", ci)})
    return out


def lane_route(spec, key_addr: str, sector):
    """The pinned route's substrate + window -- the PRODUCTION wiring,
    marked (under test): ``lane_substrate`` (the D2h landscape,
    basin-steered), ``lane_core`` (the reference's own core columns
    S-projected), ``window_columns``, the lane assembly.  Returns
    (lmol, lmf, c_core, nf, mo, ncas, act) -- ``act`` the
    PRE-projection distilled columns (the survival diagnostics read
    them; ``mo``'s window block is already projected+renormalized)."""
    from atomterm.cases import lane_substrate
    from atomterm.cells import CellKey, lane_core, window_columns
    from atomterm.windows import build_full_mo

    lmol, lmf = lane_substrate(spec)
    key = CellKey.from_address(key_addr)
    _core, nf, c_core, _rep = lane_core(key, lmol, lmf,
                                        spec.shell_bands)
    act, ncas, _w = window_columns(key, lmol, lmf, spec.shell_bands,
                                   spec.shell_radials)
    mo = build_full_mo(lmf, list(range(nf)), act, c_core=c_core)
    return lmol, lmf, c_core, nf, mo, ncas, act


def pinned_ci_run(lmol, lmf, mo, nf, ncas, sector, L, parity, c_core,
                  conv_tol: float = 1e-12, c_read=None):
    """The pinned route at the CI level (``run_irrep_casci`` -- the
    exact per-component roots at the fixed window, no optimizer).
    Returns the roots dict ``{irrep: (e_arr, ci)}`` and the columns
    the solvers read (``c_read``, the caller's symmetrized block when
    given)."""
    from atomterm.casscf import run_irrep_casci
    if c_read is None:
        c_read = np.asarray(mo)[:, nf:nf + ncas]
    return run_irrep_casci(lmol, lmf, c_read, nf, sector[0],
                           sector[1], L, parity, conv_tol=conv_tol,
                           c_core=c_core), c_read


def pinned_sa_run(lmol, lmf, mo, nf, ncas, sector, L, parity, c_core,
                  conv_tol: float = 1e-12, c_read=None,
                  mc_opts: dict = None):
    """The pinned route at the SA level (``run_irrep_sa``'s two-pass
    pattern -- the ensemble mix kernel, then the exact per-component
    roots re-solved in the frozen converged window).  Returns
    (out, c_read): the result dict (``e_tot``, ``mo``, ``roots``);
    ``mc_opts`` -- the ensemble optimizer's budget (the user's
    200-macro-cycle directive; ``None`` = PySCF's defaults)."""
    from atomterm.casscf import run_irrep_sa
    if c_read is None:
        c_read = np.asarray(mo)[:, nf:nf + ncas]
    out = run_irrep_sa(lmol, lmf, c_read, nf, sector[0], sector[1],
                       L, parity, conv_tol=conv_tol, c_core=c_core,
                       mc_opts=mc_opts)
    return out, c_read


def make_blocks(lmol, lmf, mo, nf, ncas, sector, L, parity, c_core,
                roots, c_read=None):
    """The pinned solvers' per-block kernel-ready ``CASCI`` objects on
    the given window and roots (the lane construction -- the pinned
    solver rides the block so mrpt's canonicalization stays in it;
    ``.e_tot`` carries the block's own bookkeeping).  Returns
    (blocks, expect, c_read): ``[(irrep, mc)]``."""
    from atomterm.casscf import irrep_solvers
    from atomterm.core import closed_shell_core_energy
    from atomterm.irreps import branch_d2h, orb_irreps
    from pyscf import mcscf as _mcscf

    if c_read is None:
        c_read = np.asarray(mo)[:, nf:nf + ncas]
    else:
        c_read = np.asarray(c_read)
        # the block mc must carry the columns the solvers read (the
        # symmetrized ones) in its active slots
        mo = np.array(mo, dtype=float)
        mo[:, nf:nf + ncas] = c_read
    orbsyms = orb_irreps(lmol, c_read)
    pinned = {irr: s for irr, s, _ in
              irrep_solvers(lmol, orbsyms, L, parity, conv_tol=1e-12)}
    core_part = closed_shell_core_energy(lmf, c_core)
    expect = meters(lmol, c_read, *sector)
    blocks = []
    for irrep, mult in branch_d2h(L, parity).items():
        e_arr, ci = roots[irrep]
        e_arr = np.atleast_1d(np.asarray(e_arr, dtype=float))
        mc = _mcscf.CASCI(lmf, ncas, sector)
        mc.mo_coeff = np.asarray(mo, dtype=float)
        mc.fcisolver = pinned[irrep]
        mc.orbsym = orbsyms
        mc.fcisolver.orbsym = np.asarray(orbsyms)
        ci_arr = np.asarray(ci)
        mc.ci = ([ci_arr] if ci_arr.ndim == 2
                 else [np.asarray(x) for x in ci_arr])
        mc.e_cas = e_arr[:mult]
        mc.e_tot = float(mc.e_cas[0] + core_part)
        blocks.append((irrep, mc))
    return blocks, expect, c_read


def nevpt_rows(mc, roots=None):
    """PT2 rows on a ready solution -- the engine's own canonicalization
    prelude + the eight class solves, via ``nevpt2_rows_on``."""
    from atomterm.perturb import nevpt2_rows_on
    return nevpt2_rows_on(mc, roots=roots)
