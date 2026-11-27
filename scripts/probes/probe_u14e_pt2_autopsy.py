"""The PT2 autopsy at (6,13) -- five NEVPT evaluations of ONE pinned
state-average wavefunction, printed only (no store writes).

THE QUESTION.  The pinned@pur SA rows once priced at -1580 mEh
((6,13)) where the free route prices -282 mEh.  Their readout paired
the CI vectors of the CONVERGED window with the START window's
columns.  What does each ingredient of that pairing do to E(2)?
Ranked hypotheses: (H1) the non-orthonormal orbital set inflates the
one-electron active-virtual coupling (the Sr norm's last term,
h1e_v[r,p] dm1[p,a] h1e_v[r,a], depends on the orbitals alone); (H2)
the foreign CI vector; (H3) small denominators.

THE FIVE OBJECTS, all from one optimizer run:
  (i)   as the defect built it -- the START window's active columns,
        the optimizer's virtuals, the optimizer's CI vectors;
  (ii)  the same foreign CI vectors on an ORTHONORMAL set -- the
        start window with its own complement: isolates H1 from H2;
  (iii) the same orbitals, the CI RE-SOLVED in the start window --
        the pinned route's legitimate frozen CI;
  (iv)  the optimizer's orbitals with the RUNNER's readout (the CI
        solved under the frozen-core fold) -- what the old code
        priced when it was consistent;
  (v)   the route's object: the optimizer's orbitals with their own
        readout -- the reference value.

Each object reports its birth certificate (the library's
``assert_consistent``), the active-virtual overlap, E(2) per
component and the term mean, the eight Koopmans classes (norm,
energy, the effective denominator norm/|energy|), and the Sr class
decomposed per virtual orbital into its two-body, mixed and one-body
terms (the decomposition of PySCF's own ``nevpt2.Sr``).

The library cannot build (i)-(iii) -- by design; this probe rebuilds
the defect's construction by hand to measure it.

Run from the repo root:
    OMP_NUM_THREADS=64 python scripts/probes/probe_u14e_pt2_autopsy.py \\
        | tee logs/probe_u14e_pt2_autopsy.log
"""

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "u14e"))

import numpy as np  # noqa: E402
from pyscf import fci, mcscf, mrpt  # noqa: E402
from pyscf.mrpt import nevpt2 as kernel  # noqa: E402

from atomterm.casscf import irrep_solvers, run_irrep_casci, run_irrep_sa
from atomterm.core import closed_shell_core_energy
from atomterm.irreps import orb_irreps

from vv import cases as ladder
from vv.pricing import price_classes
from vv.reference import Substrate
from vv.routes import BUDGET, State, assert_consistent
from vv.table import EH_EV

# the Koopmans classes in the kernel's order, with the kernel's call
# signature (mirrors atomterm.perturb's loop)
CLASSES = (("Sr", "Sr (-1)'", "cidms"), ("Si", "Si (+1)'", "cidms"),
           ("Sijrs", "Sijrs (0)", "eris"), ("Sijr", "Sijr (+1)", "dms"),
           ("Srsi", "Srsi (-1)", "dms"), ("Srs", "Srs (-2)", "dms"),
           ("Sij", "Sij (+2)", "dms"), ("Sir", "Sir (0)'", "dms"))


def gram_schmidt(cols, base, s) -> np.ndarray:
    """The columns S-orthonormalized on top of ``base`` (two
    refinement passes, member-preserving)."""
    cur = np.hstack(base)
    out = []
    for c in np.asarray(cols).T:
        c = np.array(c, dtype=float)
        for _ in range(2):
            c = c - cur @ (cur.T @ s @ c)
            c = c / np.sqrt(float(c @ s @ c))
        out.append(c)
        cur = np.hstack([cur, c[:, None]])
    return np.column_stack(out)


def hand_blocks(mol, mf, mo, ncore, ncas, roots, sector, term,
                conv_tol=1e-12) -> list:
    """The defect's construction, rebuilt by hand (the library refuses
    it): one kernel-ready CASCI block per D2h component on the MO
    matrix ``mo``, carrying ``roots`` -- whatever orbitals those were
    solved in."""
    orbsyms = orb_irreps(mol, mo[:, ncore:ncore + ncas])
    core_part = closed_shell_core_energy(mf, mo[:, :ncore])
    out = []
    for irrep, solver, mult in irrep_solvers(mol, orbsyms, *term, conv_tol):
        e_arr, ci = roots[irrep]
        e_arr = np.atleast_1d(np.asarray(e_arr, dtype=float))
        mc = mcscf.CASCI(mf, ncas, sector)
        mc.mo_coeff = mo
        mc.fcisolver = solver
        mc.orbsym = orbsyms
        mc.fcisolver.orbsym = np.asarray(orbsyms)
        ci_arr = np.asarray(ci)
        mc.ci = (ci_arr[0] if ci_arr.ndim == 3 else ci_arr) if mult == 1 \
            else [np.asarray(x) for x in ci_arr]
        mc.e_cas = e_arr[:mult]
        mc.e_tot = float(mc.e_cas[0] + core_part)
        out += [State(mc, k, irrep) for k in range(mult)]
    return out


def price(state) -> tuple:
    """(the eight classes as (label, norm, energy), the Sr
    decomposition) -- the prelude of atomterm.perturb._root_row
    (canonicalization, the 3-RDM, the two 4-pdm contractions, the
    ERIs), then the kernel's classes with Sr opened up."""
    mc, root = state.mc, state.root
    pt = mrpt.NEVPT(mc, root=root, density_fit=False)
    if pt.fcisolver.nroots > 1:
        pt.ci = list(pt.ci)
    pt.mo_coeff, ci_vec, pt.mo_energy = pt.canonicalize(
        pt.mo_coeff, ci=pt.load_ci(), cas_natorb=True, verbose=0)
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
    eris = kernel._ERIS(pt, pt.mo_coeff)
    aaaa = eris["ppaa"][ncore:nocc, ncore:nocc].copy()
    dms["f3ca"] = kernel._contract4pdm("NEVPTkern_cedf_aedf", aaaa, ci,
                                       ncas, pt.nelecas, (link_a, link_b))
    dms["f3ac"] = kernel._contract4pdm("NEVPTkern_aedf_ecdf", aaaa, ci,
                                       ncas, pt.nelecas, (link_a, link_b))
    classes = []
    for fname, label, kind in CLASSES:
        fn = getattr(kernel, fname)
        if kind == "cidms":
            norm, energy = fn(pt, ci, dms, eris)
        elif kind == "dms":
            norm, energy = fn(pt, dms, eris)
        else:
            norm, energy = fn(pt, eris)
        classes.append((label, float(norm), float(energy)))
    return classes, sr_parts(pt, dms, eris)


def sr_parts(pt, dms, eris) -> dict:
    """The Sr norm and energy per virtual orbital, each in its three
    terms (PySCF's ``nevpt2.Sr``, the eris branch, opened up)."""
    ncore, ncas = pt.ncore, pt.ncas
    nocc = ncore + ncas
    dm1, dm2, dm3 = dms["1"], dms["2"], dms["3"]
    h1e = eris["h1eff"][ncore:nocc, ncore:nocc]
    h2e = eris["ppaa"][ncore:nocc, ncore:nocc].transpose(0, 2, 1, 3)
    h2e_v = eris["ppaa"][nocc:, ncore:nocc].transpose(0, 2, 1, 3)
    h1e_v = eris["h1eff"][nocc:, ncore:nocc] \
        - np.einsum("mbbn->mn", h2e_v)
    ci = pt.load_ci()
    a16 = kernel.make_a16(h1e, h2e, dms, ci, ncas, pt.nelecas)
    a17 = kernel.make_a17(h1e, h2e, dm2, dm3)
    a19 = kernel.make_a19(h1e, h2e, dm1, dm2)
    ener = (np.einsum("ipqr,pqrabc,iabc->i", h2e_v, a16, h2e_v)
            + np.einsum("ipqr,pqra,ia->i", h2e_v, a17, h1e_v) * 2.0
            + np.einsum("ip,pa,ia->i", h1e_v, a19, h1e_v))
    two = np.einsum("ipqr,rpqbac,iabc->i", h2e_v, dm3, h2e_v)
    mixed = np.einsum("ipqr,rpqa,ia->i", h2e_v, dm2, h1e_v) * 2.0
    one = np.einsum("ip,pa,ia->i", h1e_v, dm1, h1e_v)
    return {"two": two, "mixed": mixed, "one": one, "ener": ener,
            "diff": np.asarray(pt.mo_energy[nocc:]),
            "h1e_v2": np.einsum("ip,ip->i", h1e_v, h1e_v)}


def report(name, states, s, ncore, ncas) -> dict:
    print(f"\n===== {name} =====", flush=True)
    mo = states[0].mo_coeff
    act, virt = mo[:, ncore:ncore + ncas], mo[:, ncore + ncas:]
    overlap = np.abs(act.T @ s @ virt)
    print(f"  active-virtual overlap max|A^T S V| {overlap.max():.3e}",
          flush=True)
    try:
        cert = assert_consistent(states, None, "averaged")
        print(f"  certificate: PASSES -- orbitals {cert['ortho']:.1e}, "
              f"eigenpairs {cert['residual']:.1e}", flush=True)
    except ValueError as e:
        print(f"  certificate: REFUSED -- {e}", flush=True)
    t0 = time.time()
    e2, cls, sr = [], {}, []
    for st in states:
        classes, parts = price(st)
        e2.append(sum(e for _l, _n, e in classes))
        for label, norm, energy in classes:
            n, e = cls.get(label, (0.0, 0.0))
            cls[label] = (n + norm, e + energy)
        sr.append(parts)
    e2 = np.array(e2)
    print(f"  E(2) per component (mEh): "
          + " ".join(f"{x * 1e3:+.1f}" for x in e2), flush=True)
    print(f"  E(2) term mean {e2.mean() * 1e3:+.4f} mEh "
          f"({e2.mean() * EH_EV:+.4f} eV); spread {np.ptp(e2) * 1e3:.2f} mEh "
          f"({time.time() - t0:.0f} s)", flush=True)
    print("  class          norm(mean)   energy(mean, mEh)   eff. denom.", flush=True)
    for _fname, label, _kind in CLASSES:
        norm, energy = (cls[label][0] / len(states), cls[label][1] / len(states))
        eff = norm / abs(energy) if energy else float("nan")
        print(f"  {label:12s} {norm:11.4f}   {energy * 1e3:+12.4f}      "
              f"{eff:9.3f}", flush=True)
    one = np.mean([p["one"].sum() for p in sr])
    mixed = np.mean([p["mixed"].sum() for p in sr])
    two = np.mean([p["two"].sum() for p in sr])
    print(f"  Sr norm decomposed (mean over components): one-body {one:.4f} "
          f"+ mixed {mixed:.4f} + two-body {two:.4f}", flush=True)
    print(f"  max|h1e_v[r,a]|^2 per virtual (mean): "
          f"{np.mean([p['h1e_v2'].max() for p in sr]):.3e}", flush=True)
    worst = np.argsort(-np.mean([p["one"] for p in sr], axis=0))[:4]
    for r in worst:
        print(f"    virtual {r:3d}: one-body norm "
              f"{np.mean([p['one'][r] for p in sr]):9.4f}, denominator "
              f"{np.mean([p['diff'][r] for p in sr]):8.4f}", flush=True)
    return {"e2": e2.mean(), "sr_one": one, "overlap": overlap.max()}


def main() -> int:
    t00 = time.time()
    sub = Substrate.load()
    built = ladder.CASES["(6,13)"].window(sub)
    space = built.space
    space.certify(sub.mf)
    for line in built.report:
        print(f"  {line}", flush=True)
    mol, mf = sub.d2h
    s = mf.get_ovlp()
    ncore, ncas, sector, term = space.ncore, space.ncas, space.sector, ladder.TERM
    na, nb = sector

    print("\n===== the optimizer run (the pinned SA, the production runner) =====",
          flush=True)
    t0 = time.time()
    out = run_irrep_sa(mol, mf, space.active, ncore, na, nb, *term,
                       conv_tol=1e-12, c_core=space.core, mc_opts=BUDGET)
    print(f"  ensemble {out['e_tot']:.8f}, converged {out['converged']} "
          f"(cycles {out['cycles']}, conv dE {out['conv_dE']:.1e}, "
          f"{time.time() - t0:.0f} s)", flush=True)
    mo_c = np.asarray(out["mo"])
    start, conv = space.active, mo_c[:, ncore:ncore + ncas]
    cos = np.linalg.svd(start.T @ s @ conv, compute_uv=False)
    print(f"  the start window vs the converged active block: cosines "
          f"[min {cos.min():.4f} / max {cos.max():.4f}], shared directions "
          f"({np.sum(cos > 0.99)}/{ncas})", flush=True)

    # the orbitals of the orthonormal objects: the start window with
    # its own complement
    ortho = np.hstack([mo_c[:, :ncore], start,
                       gram_schmidt(mo_c[:, ncore + ncas:],
                                    [mo_c[:, :ncore], start], s)])
    as_bugged = np.array(mo_c)
    as_bugged[:, ncore:ncore + ncas] = start
    print("\n===== the readouts =====", flush=True)
    t0 = time.time()
    roots_frozen = run_irrep_casci(mol, mf, start, ncore, na, nb, *term,
                                   conv_tol=1e-12, c_core=space.core)
    print(f"  the CI re-solved in the START window ({time.time() - t0:.0f} s)",
          flush=True)
    roots_own = run_irrep_casci(mol, mf, conv, ncore, na, nb, *term,
                                conv_tol=1e-12, c_core=mo_c[:, :ncore])

    objects = (
        ("(i) as the defect built it: start window, optimizer virtuals, "
         "optimizer CI", as_bugged, out["roots"]),
        ("(ii) orthonormal: start window + its complement, the SAME foreign CI",
         ortho, out["roots"]),
        ("(iii) orthonormal: start window + its complement, the CI re-solved "
         "in it", ortho, roots_frozen),
        ("(iv) the optimizer's orbitals, the RUNNER's readout (frozen-core "
         "fold)", mo_c, out["roots"]),
        ("(v) CONSISTENT: the optimizer's orbitals, their own readout",
         mo_c, roots_own),
    )
    results = {}
    for name, mo, roots in objects:
        states = hand_blocks(mol, mf, mo, ncore, ncas, roots, sector, term)
        results[name] = report(name, states, s, ncore, ncas)

    # the instrument's self-check: the class sum against the engine
    # on the consistent object
    states = hand_blocks(mol, mf, mo_c, ncore, ncas, roots_own, sector, term)
    mine = sum(e for _l, _n, e in price(states[0])[0])
    engine = price_classes(states[0]).e_corr
    print(f"\n===== the instrument =====\n  the class sum against the "
          f"engine's PT2Row on component 0 of (v): {mine:.12f} vs "
          f"{engine:.12f} Eh (diff {mine - engine:.1e})", flush=True)

    print("\n===== the verdict of the measurements =====", flush=True)
    ref = results[objects[-1][0]]["e2"]
    for name, r in results.items():
        print(f"  {name.split(':')[0]:52s} E(2) {r['e2'] * 1e3:+10.2f} mEh "
              f"(vs consistent {abs(r['e2'] - ref) * 1e3:9.2f} mEh off); "
              f"Sr one-body {r['sr_one']:8.4f}; overlap {r['overlap']:.3f}",
              flush=True)
    print(f"\n({time.time() - t00:.0f} s total)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
