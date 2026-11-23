"""The u14e V&V family's shared core: metadata, the trusted rungs, the
committed anchors, the analysis meters, the derivations' thin
helpers, the verdict ledger, the row store and the one rendered
table.  The DERIVATIONS live in the legs -- each leg readable
top-to-bottom; the free routes are PySCF-direct (independent of the
production construction chain), the pinned routes and the 6_13
distilled radial are marked production calls (under test /
genuinely distilled machinery).

Run the legs from the repo root (or all through the driver):
    python scripts/u14e/u14e_vv_4_12.py [--ground-only]
    python scripts/u14e/u14e_vv_4_15.py [--ground-only]
    python scripts/u14e/u14e_vv_6_13.py [--ground-only]
    python scripts/u14e/u14e_vv_uhf.py
    python scripts/u14e/u14e_vv_orchestration.py [--ground-only]
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

import numpy as np  # noqa: E402
from pyscf import mcscf, mrpt, scf  # noqa: E402

from atomterm.cases import campaign_spec  # noqa: E402
from atomterm.detbasis import determinant_basis, operator_apply  # noqa: E402
from atomterm.operators import (  # noqa: E402
    ci_in_detbasis_order, l2_operator_parts, s2_operator_parts)

EH_EV = 27.211386245988
L_L2 = 72.0  # the 5L parent term: L = 8 -> L(L+1)
S_S2 = 6.0   # spin 4 -> S = 2 -> S(S+1)

ROWS = Path(__file__).resolve().parent / "vv.json"
TABLE = Path(__file__).resolve().parent / "VV.md"
PAYLOAD = Path("elements") / "U14e" / "perturbation.json"
RUNGS = Path("elements") / "U14e" / "rungs"

# the assertion bands (declared constants; each cites the observed
# fresh values -- logs/u14e_vv_*.log):
BAND_L2 = 1e-5     # the Casimir meters are gauge-invariant; the
                   # seeded readout's observed dev 3.0e-06 (free)
                   # and 1.9e-07 (pinned, all 72.000000)
BAND_E2_MEH = 2.0  # a fresh SA draw is a new orbital gauge and E(2)
                   # is state-specific; observed delta 0.0006-0.0029
                   # mEh (the engine wobble is ~0.3 uEh run-to-run)
BAND_ROUTE = 1e-6  # Eh; the ENSEMBLE energy is the route comparison
                   # surface (variational, gauge-invariant; E(2) is
                   # orbital-sensitive and printed, never asserted);
                   # observed 3.3e-09 Eh
BAND_CLASS = 0.01  # the certified-class meters (the (6,13) pinned
                   # components: 16 at L(L+1)=90 + one at 72)


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


# the analysis meters --------------------------------------------------

def meters(mol, c_act, na: int, nb: int):
    """The analysis meters as their ACTION (atomterm's angular
    operators through the matrix-free application; the many-body
    meters read per root through the CI-order bridge).  The one piece
    of atomterm analysis the legs keep -- the construction is
    double-certified in angular.py; the meters are the terms'
    identity, not part of the derivation under test."""
    norb = c_act.shape[1]
    dets = determinant_basis(norb, na, nb)
    apply_l2 = operator_apply(dets, *l2_operator_parts(mol, c_act))
    apply_s2 = operator_apply(dets, *s2_operator_parts(norb))
    apply = {"l2": apply_l2, "s2": apply_s2}

    def expect(kind: str, ci_pyscf) -> float:
        w = ci_in_detbasis_order(np.asarray(ci_pyscf).ravel(),
                                 norb, na, nb)
        return float(np.real(w @ apply[kind](w)))

    return expect


# the derivations' thin helpers (one home; the legs show the setup) ----

def free_kernel(mf, ncas: int, sector, mo_full, weights=None,
                nroots: int = 1, track_iters: bool = False):
    """The free route's two kernels, written once: CASSCF (plain
    single root, or ``state_average_`` over ``weights``), then the
    two-pass readout (CASCI, ``canonicalization = False`` -- the
    window exactly as converged), seeded from the SA kernel's own
    vectors (the faithful two-pass form; unseeded, the davidson tails
    mix a quasi-degenerate block at the 1e-5 level and the meters
    carry it).  Returns (mc, cas).  ``track_iters`` rides the
    optimizer's callback (its ``locals()`` carry ``imacro``/``de``/
    ``conv``) and sets ``mc.n_macro``/``mc.final_de`` -- the conv
    record; the default off leaves the historical path untouched."""
    mc = mcscf.CASSCF(mf, ncas, sector)
    if track_iters:
        _it = {"imacro": -1, "de": float("nan"), "conv": None}

        def _cb(envs):
            if envs.get("imacro", -1) > _it["imacro"]:
                _it["imacro"] = envs["imacro"]
                _it["de"] = abs(float(envs.get("de", float("nan"))))
                _it["conv"] = envs.get("conv")
        mc.callback = _cb
    if weights is not None and len(weights) > 1:
        mc.state_average_(list(weights))
    mc.kernel(mo_full)
    if track_iters:
        mc.n_macro = _it["imacro"] + 1
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
                  conv_tol: float = 1e-12, c_read=None):
    """The pinned route at the SA level (``run_irrep_sa``'s two-pass
    pattern -- the ensemble mix kernel, then the exact per-component
    roots re-solved in the frozen converged window).  Returns
    (out, c_read): the result dict (``e_tot``, ``mo``, ``roots``)."""
    from atomterm.casscf import run_irrep_sa
    if c_read is None:
        c_read = np.asarray(mo)[:, nf:nf + ncas]
    out = run_irrep_sa(lmol, lmf, c_read, nf, sector[0], sector[1],
                       L, parity, conv_tol=conv_tol, c_core=c_core)
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


def _irrep_projectors(mol, S):
    """``[(D_k, G_k^{-1} D_k^T S)]`` per irrep block -- the metric
    projector pieces, ``P_k v = D_k @ (Ginv_dTS @ v)``.  The irrep
    spans are mutually S-orthogonal and partition the AO space, so
    ``sum_k ||P_k v||^2_S = ||v||^2_S`` whatever the blocks' internal
    gauge (pyscf keeps raw AO selectors in single-irrep blocks,
    max|D^T S D - I| = 0.25 on the sto-3g probe)."""
    projs = []
    for d in mol.symm_orb:
        d = np.asarray(d)
        G = d.T @ S @ d
        projs.append((d, np.linalg.solve(G, d.T @ S)))
    return projs


def symmetrize_columns(mol, S, c, tol: float = 1e-4):
    """Re-project each column onto its DOMINANT D2h irrep -- the
    disclosed CLEANUP for the tiny irrep tails a core projection's
    conditioning leaves (``orb_irreps``/``label_orb_symm`` refuses
    past ~1e-6); the core projection itself is physics and is never
    undone here.  The projection is the METRIC projector
    ``P_k = D_k (D_k^T S D_k)^{-1} D_k^T S`` per irrep block: the
    irrep spans are mutually S-orthogonal and partition the AO space,
    so ``sum_k ||P_k v||^2_S = ||v||^2_S`` exactly and the impurity
    ``1 - max_k ||P_k v||^2_S / ||v||^2_S`` sits in ``[0, 1]``
    whatever the blocks' internal gauge (a coordinate-only
    projection overcounts when ``D_k`` is not S-orthonormal -- the
    defect the certification caught).  The impurity is measured on
    the S-normalized column, so it is scale-free; already-pure
    columns project to themselves (the identity on their span) and
    come back just normalized.  Returns ``(c_fixed, impurities)`` --
    the impurities are the per-column BEFORE values, the caller
    discloses them; raises past ``tol`` (a real mixing class, not
    the projection-conditioning dirt this cleanup owns)."""
    projs = _irrep_projectors(mol, S)
    fixed = np.empty_like(c)
    imp = []
    for j in range(c.shape[1]):
        v = np.asarray(c[:, j], dtype=float)
        n2 = float(v @ S @ v)
        best_w, best_n = None, -1.0
        for d, Ginv_dT in projs:
            w = d @ (Ginv_dT @ v)
            n_k = float(w @ S @ w)
            if n_k > best_n:
                best_w, best_n = w, n_k
        imp.append(1.0 - best_n / n2)
        fixed[:, j] = best_w / np.sqrt(best_n)
    worst = max(imp)
    if worst > tol:
        raise ValueError(
            f"a column carries irrep impurity {worst:.2e} past the "
            f"tolerance {tol:.1e} -- a real mixing class, not the "
            f"projection-conditioning dirt this cleanup owns")
    return fixed, imp


def dominant_irreps(mol, S, c):
    """Per column, the id of the irrep span holding the most of it
    (the metric-projector weights); :func:`symmetrize_columns`
    re-projects onto these."""
    projs = _irrep_projectors(mol, S)
    ids = []
    for j in range(c.shape[1]):
        v = np.asarray(c[:, j], dtype=float)
        best_k, best_n = 0, -1.0
        for k, (d, Ginv_dTS) in enumerate(projs):
            w = d @ (Ginv_dTS @ v)
            nk = float(w @ S @ w)
            if nk > best_n:
                best_k, best_n = k, nk
        ids.append(best_k)
    return ids


def purify_window(mol, S, c_core, c, tol: float = 1e-10):
    """Converge a window block to the INTERSECTION of the frozen-core
    S-orthogonalization and the dominant-irrep purity, in closed form.

    The naive iterate (project core, project irrep, renormalize)
    converges only geometrically on a core-lean column (measured at
    (4,15): rate ~0.36 at 92% core-lean -- stalled at 5e-04 after 8
    passes; a single pass leaves a 2.7e-07 tail past pyscf's orb
    labeling).  The closed form removes the core's content
    within each dominant-irrep block exactly: the block subtracts
    ``B_k = P_k(c_core)`` -- ALL core columns projected into
    ``span(D_k)`` -- so no assignment of core columns to blocks is
    needed (the lane core's global orthonormalization mixes across
    irrep blocks; a per-assigned-block form stalls at 5.3e-04).

    Returns ``(c_fixed, report)`` with ``report = (passes,
    core_res, irrep_imp, ortho_dev)``; raises on genuinely-mixed
    input (via :func:`symmetrize_columns`' contract), on a vanished
    column (the column lies INSIDE the frozen core's span -- it has
    no core-orthogonal remainder at all; the measured (4,15) class:
    the C1 reference's m-mixed occupied core, S-projected onto the
    declared shell span, contributes the promoted band's own slots),
    or when the verification residues miss ``tol``.  THE CONTRACT:
    the caller passes the ASSEMBLED block (the pipeline's form --
    ``build_full_mo``'s output, already core-orthogonal); a
    genuinely-core-lean RAW column is refused -- core-projecting it
    is the assembly's job, and a global pre-projection here would
    SYNTHESIZE irrep content against a cross-irrep-mixed core (the
    defect the certification caught: the remainder picks up the
    core's other-block directions at the correlation scale).  THE
    CORE NEED NOT BE IRREP-PURE: the lane core's global
    orthonormalization mixes across irrep blocks (measured:
    per-assigned-block subtraction stalls at a 5.3e-04 core
    residue), so each block removes the core's COMPLETE content
    within its span, ``B_k = P_k(c_core)`` -- assignments
    irrelevant."""
    v, imp = symmetrize_columns(mol, S, c)
    ids = dominant_irreps(mol, S, v)
    projs = _irrep_projectors(mol, S)
    fixed = np.empty_like(v)
    for k in sorted(set(ids)):
        cols = [j for j, kid in enumerate(ids) if kid == k]
        d, Ginv_dTS = projs[k]
        sub = v[:, cols]
        B = d @ (Ginv_dTS @ c_core)     # ALL core columns, in span(D_k)
        Gb = B.T @ S @ B
        gb, Ub = np.linalg.eigh(Gb)
        keep = gb > 1e-10 * gb.max()
        Bk = B @ Ub[:, keep] / np.sqrt(gb[keep])
        sub = sub - Bk @ (Bk.T @ S @ sub)   # the block's complete core
        sub = d @ (Ginv_dTS @ sub)          # the irrep cleanup: no-op
        norms = np.sqrt(np.diag(sub.T @ S @ sub))
        if norms.min() < 1e-8:
            raise VanishedColumnError(
                f"a window column vanished into the core (post-"
                f"subtraction S-norm {norms.min():.1e}) -- the column "
                f"lies inside the frozen core's span and has no core-"
                f"orthogonal remainder")
        sub = sub / norms
        g = sub.T @ S @ sub
        gv, U = np.linalg.eigh(g)
        if gv.min() <= 1e-10 * max(gv.max(), 1.0):
            raise ValueError(
                f"the window block is rank-deficient in irrep {k} "
                f"(smallest Gram eigenvalue {gv.min():.1e}) -- "
                f"dependent columns, not a window")
        sub = sub @ U @ np.diag(gv ** -0.5) @ U.T
        for q, j in enumerate(cols):
            fixed[:, j] = sub[:, q]
    core_res = float(np.linalg.norm(c_core.T @ S @ fixed, axis=0).max())
    _, imp_v = symmetrize_columns(mol, S, fixed, tol=1.0)
    ortho_dev = float(np.max(np.abs(
        fixed.T @ S @ fixed - np.eye(fixed.shape[1]))))
    if core_res > tol or max(imp_v) > tol or ortho_dev > 1e-8:
        raise ValueError(
            f"the closed-form intersection missed its tolerances "
            f"(core residue {core_res:.1e}, irrep impurity "
            f"{max(imp_v):.1e}, orthonormality {ortho_dev:.1e})")
    return fixed, {"passes": 1, "core_res": core_res,
                   "irrep_imp": float(max(imp_v)), "ortho_dev": ortho_dev}


def meters_all(blocks, expect):
    """The meters over EVERY block component (in branch order) --
    independent of any PT2 root cut."""
    l2s, s2s = [], []
    for _, mc in blocks:
        for k in range(len(np.atleast_1d(mc.e_cas))):
            ci = np.asarray(mc.ci[k]).ravel()
            l2s.append(expect("l2", ci))
            s2s.append(expect("s2", ci))
    return l2s, s2s


# the verdict ledger ---------------------------------------------------


class VanishedColumnError(ValueError):
    """A window column vanished into the core inside
    :func:`purify_window` -- the column lies inside the frozen core's
    span (no core-orthogonal remainder exists for it; the measured
    (4,15) class: the C1 reference's m-mixed occupied core, projected
    onto the declared shell span, spans the promoted band's own
    slots)."""

class Verdict:
    """Every check lands one row; the exit code is the verdict.
    ``note`` rows are measured and recorded, never gating (the
    new-physics class)."""

    def __init__(self):
        self.rows = []

    def check(self, name, ok, detail):
        ok = bool(ok)
        self.rows.append((name, ok, "check"))
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}",
              flush=True)
        return ok

    def note(self, name, detail):
        self.rows.append((name, True, "note"))
        print(f"  [note] {name}: {detail}", flush=True)

    def report(self) -> bool:
        print("\n===== verdict =====", flush=True)
        for name, ok, kind in self.rows:
            tag = ("PASS" if ok else "FAIL") if kind == "check" \
                else "REPORT"
            print(f"  {tag}  {name}", flush=True)
        checks = [(n, ok) for n, ok, k in self.rows if k == "check"]
        good = sum(ok for _, ok in checks)
        print(f"{good}/{len(checks)} assertions hold", flush=True)
        return good == len(checks)


# the row store + the one table ---------------------------------------

def de_fmt(de_eh: float) -> str:
    """One delta, both units: ``mEh / eV`` (4 decimals)."""
    return f"{de_eh * 1e3:+.4f} / {de_eh * EH_EV:+.4f}"


def drop_row(row: dict) -> None:
    """One leg's result row, merged into the family's store
    (``method`` + ``route`` + ``cas`` is the key)."""
    rows = json.loads(ROWS.read_text()) if ROWS.exists() else {}
    rows[f"{row['method']}|{row['route']}|{row['cas']}"] = row
    ROWS.write_text(json.dumps(rows, indent=1) + "\n")
    render_table()


_ORDER = (["rohf|-|-", "uhf|-|-"] +
          [f"{m}|{r}|{c}"
           for m, r in (("casci", "free"), ("casscf", "free"),
                        ("sa-casscf", "free"), ("casci", "pinned"),
                        ("sa-casscf", "pinned"))
           for c in ("(4,12)", "(4,15)", "(6,13)")])


def fmt_range(vals, nd: int = 4) -> str:
    """One value when the set agrees, ``a..b`` when it scatters."""
    lo, hi = min(vals), max(vals)
    if hi - lo < 10 ** (-nd) / 2:
        return f"{lo:.{nd}f}"
    return f"{lo:.{nd}f}..{hi:.{nd}f}"


def conv_of(precision, n_iters, converged=True) -> str:
    """The conv cell: ``"<precision>/<n_iters>"`` (the user's form).
    ``precision`` -- the achieved |dE| where the optimizer reports it
    (PySCF's callback ``de``; ``run_irrep_sa``'s ``conv_dE``), else
    the requested tolerance; ``n_iters`` -- the optimizer's macro
    cycle count, ``None`` for the CI solves (no iteration
    structure); ``*`` marks a cycle-cap draw (not converged)."""
    star = "" if converged else "*"
    p = "—" if precision is None else f"{precision:.1e}"
    n = "—" if n_iters is None else str(n_iters)
    return f"{p}/{n}{star}"


def render_table() -> None:
    """The small committed table: the reference rows lead (ROHF is
    the zero of every dE), then one row per method x route x CAS.
    Every energy <= 4 decimals; deltas in the ``mEh / eV`` slash
    form."""
    rows = json.loads(ROWS.read_text()) if ROWS.exists() else {}
    order = [k for k in _ORDER if k in rows]
    order += [k for k in rows if k not in order]
    e0 = rohf_energy()
    out = [
        "# U14e -- the V&V family",
        "",
        "| method | route | CAS | &lt;L²&gt; | &lt;S²&gt; | conv "
        "| E_method (Eh) | dE_method (mEh/eV) | E(2) (mEh/eV) "
        "| E_tot (Eh) | dE_tot (mEh/eV) | term |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for k in order:
        r = rows[k]
        de_var = r["e_method"] - e0
        de_tot = r["e_tot"] - e0
        e2 = (f"{r['e2_meh']:+.4f} / {r['e2_ev']:+.4f}"
              if r.get("e2_meh") is not None else "—")
        # the render-time consistency check: the store's own bookkeeping
        if r.get("e2_meh") is not None and abs(
                r["e_method"] + r["e2_meh"] * 1e-3 - r["e_tot"]) > 5e-4:
            raise ValueError(
                f"render check: row {k} violates E_tot = E_method + "
                f"E(2) by {r['e_method'] + r['e2_meh'] * 1e-3 - r['e_tot']:.2e} Eh")
        out.append(
            f"| {r['method']} | {r['route']} | {r['cas']} "
            f"| {r['l2']} | {r['s2']} | {r.get('conv') or '—'} "
            f"| {r['e_method']:.4f} | {de_fmt(de_var)} | {e2} "
            f"| {r['e_tot']:.4f} | {de_fmt(de_tot)} "
            f"| {r.get('term', '—')} |")
    out += [
        "",
        "E_method -- the method's own energy (the CASSCF/SA total; "
        "the pinned casci value is the TERM energy, the components' "
        "mean plus the core bookkeeping constant -- its dE is that "
        "constant's projected-density form, not comparable to the "
        "free rows).  E_tot = E_method + E(2) (SC-NEVPT2, pyscf.mrpt, "
        "strongly contracted); dE against the ROHF row.",
        "",
        "E(2) comparability -- within a row, E_tot = E_method + E(2); "
        "across rows, compare E_tot at matched partitions.  The "
        "frozen-window rows' large E(2) is the UNOPTIMIZED gauge: "
        "the optimizer converts PT2-priced (dynamic) correlation "
        "into variationally-absorbed (static) correlation -- read "
        "that trade down the E_method/E(2) columns (at (4,12): "
        "casci E(2) -1054 mEh -&gt; SA -342 mEh while E_var falls "
        "+923 -&gt; -42 mEh; the net is the E_tot column).",
        "",
        "&lt;L²&gt;/&lt;S²&gt; -- the ground root for the single-root "
        "methods; min..max over the ensemble's 17 roots (free SA) or "
        "over the term's pinned components (pinned rows).  The spread "
        "is the disclosure (draw-dependent partitions, closure costs, "
        "block-ground classes); the per-root record lives in "
        "`logs/u14e_vv_*.log`, the committed campaign payload in "
        "`elements/U14e/perturbation.json` -- two independent "
        "sources, compared here.",
        "",
        "conv -- \"&lt;precision&gt;/&lt;n_iters&gt;\": the "
        "optimizer's achieved |dE| (the callback's final macro-step "
        "delta; the pinned SA's conv_dE) over its macro-cycle count; "
        "the CI solves carry the requested CI tolerance with no "
        "iteration structure (—).  ``*`` marks a cycle-cap draw "
        "(not converged; the run's log discloses).  Cells read — "
        "where a run's record predates the column (the convergence "
        "flag lives in that run's log).",
        "",
    ]
    if "casci|pinned|(4,15)" in rows:
        out += ["The pinned (4,15) row is the RE-ATTEMPT under the "
                "corrected symmetrizer, and it MEASURED the refusal: "
                "the declared core span of @6s6p7s realizes 6p as "
                "THREE bands (13p/14p/15p), and the pinned route's "
                "frozen core (the C1 reference's m-mixed occupied 6p "
                "projected onto that span) spans the 14p slots "
                "ENTIRELY -- the promoted columns are core-span "
                "content by construction, the pinned SA has no "
                "window at this partition (the typed "
                "VanishedColumnError, re-measured every run), and "
                "the pinned casci's block grounds come out NOT the "
                "5L (L classes 0..42; at (4,12) the same route reads "
                "17 components at 72.000000).  The pinned (4,15) is "
                "the pinned (4,12) here.  The FREE route's core (the "
                "pur substrate) does not contain the slots (survival "
                "0.28) and hosts the promotion -- its measured cost "
                "+3.09 mEh (fresh, twice) is the experiment's "
                "answer; the lane_core promotion is declined, src "
                "unchanged.", ""]
    else:
        out += ["The pinned (4,15) legs are measured out: the "
                "declared core span of @6s6p7s realizes 6p as three "
                "bands (13p/14p/15p), and the pinned route's frozen "
                "core spans the 14p slots entirely -- the promoted "
                "columns are core-span content by construction (the "
                "typed VanishedColumnError).  The free route hosts "
                "the promotion.", ""]
    if "casci|free|(6,13)" in rows:
        out += ["The (6,13) rows carry the CORRECTED meters (the "
                "2026-11-22 slice incident: the leg's meter slices "
                "ran one column off the key's core-pair count; the "
                "kernels were never affected).  The frozen ground is "
                "the 7F -- L = 3, S = 3 definite (the d5 high-spin "
                "class + one f; the 7s column UNOCCUPIED at 3e-03): "
                "the as-built window's variational ground abandons "
                "the intended 5f3 6d1 7s2 configuration, driven by "
                "the weak distillation captures (10f 0.82, 13d 0.92, "
                "7s 0.27-0.71).  The window-rework candidate is the "
                "s distillation's k (2 or 3 of the 3 s bands) -- a "
                "declaration change, the user's call.  The window "
                "itself is healthy: assembled span = projected span "
                "exactly, closure meters ~1e-13 (the campaign cell "
                "agreed all along).", ""]
    TABLE.write_text("\n".join(out))


def std_row(method: str, route: str, cas: str, e_method: float,
            e2_meh: float | None, e_tot: float, l2: str, s2: str,
            term: str = "—", conv: str | None = None) -> dict:
    """One table row in the family's schema (the deltas are computed
    at render time against the ROHF row; ``l2``/``s2`` pre-formatted;
    ``e_method`` -- the method's own energy, whatever its semantics;
    ``e_tot`` -- the PT2-corrected total, ``e_method + E(2)``;
    ``term`` -- the row's term-identity disclosure, the pick
    column; ``conv`` -- the ``"<precision>/<n_iters>"`` cell,
    ``conv_of``'s form)."""
    return {
        "method": method, "route": route, "cas": cas,
        "l2": l2, "s2": s2,
        "e_method": e_method,
        "e2_meh": e2_meh,
        "e2_ev": None if e2_meh is None else e2_meh * EH_EV / 1e3,
        "e_tot": e_tot,
        "term": term,
        "conv": conv,
    }
