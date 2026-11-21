"""The CAS(4,12) V&V leg -- the 5L state-averaged leg, PySCF-direct
and readable top-to-bottom.

THE SETUP (all visible below):
  1. the trusted ROHF rung LOADED (never re-run -- the cache is the
     reproducibility unit); its identity cross-checked against the
     campaign declaration;
  2. the core: the rotation-closed span over the occupied prefix
     (atomterm.core's certified construction -- a raw core polarizes
     and splits the term, and the committed anchors live on this
     substrate; this is substrate chemistry, not the derivation
     under test);
  3. the window: the AO columns of the 10f (7) + 13d (5) shells,
     written out -- selected by label prefix, projected out of the
     core span, Loewdin-orthonormalized, the S-orthonormal
     complement filled by the deflated overlap.  These shells are
     band-pure at k=1, so the visible selection IS the stored
     distilled window; the production builder's product is
     cross-checked elementwise (the certification assert);
  4. the SA: PySCF's CASSCF with state_average_ over the ground
     term's 17 M_L components (L = 8, u);
  5. the two-pass readout (CASCI, canonicalization off) in the
     frozen converged window;
  6. PT2: the engine's own mrpt.NEVPT(mc, root=k).kernel() per root;
  7. the meters: atomterm's exact angular operators (shared core).

THE ASSERTIONS:
  (a) fresh E(2) root 0 ~ the committed payload row
      (u14e_10f13d_4-0_sa-casscf_ao_-_pur_avg, root 0), band +-2 mEh;
  (b) <L2> = 72 / <S2> = 6 on every root (band 1e-5);
  --pinned adds (d): the D2h-pinned route ~ this visible free route.
  The pinned route is the PRODUCTION evaluator (run_irrep_sa on the
  lane substrate) -- it is the thing under test; the free route
  above is the independent reference.

Run from the repo root:
    OMP_NUM_THREADS=16 python scripts/u14e/u14e_vv_4_12.py \
        [--pinned] [--ground-only]  |  tee logs/u14e_vv_4_12.log
"""

import sys
import time
from datetime import date

import numpy as np
from pyscf import gto, mcscf, mrpt, scf

from _vv_common import (ATOM_DIR, BAND_E2_MEH, BAND_L2, BAND_ROUTE,
                        EH_EV, L_L2, S_S2, Verdict, anchor, drop_row,
                        load_reference, meters)

ADDR_412 = "u14e_10f13d_4-0_sa-casscf_ao_-_pur_avg"
SECTOR = (4, 0)   # (na, nb): the maximal-Ms sector, 4 open alpha
NCORE = 5         # 6s2 6p6 7s2 -- five closed pairs
SHELLS = (("10f", 7), ("13d", 5))   # the window's AO shells


def visible_window(mol, mf, c_core):
    """The window, written out: select the shells' AO columns, project
    the core span out (the S metric), Loewdin-orthonormalize, fill the
    complement by the deflated overlap.  Returns (mo_full, c_act) and
    asserts elementwise agreement with the production builder."""
    labels = [lab.split()[-1] for lab in mol.ao_labels()]
    sel = [i for i, lab in enumerate(labels)
           if any(lab.startswith(w) for w, _ in SHELLS)]
    ncas = len(sel)
    if ncas != sum(n for _, n in SHELLS):
        sys.exit(f"the label selection found {ncas} AOs, expected "
                 f"{sum(n for _, n in SHELLS)}")
    nao = mol.nao
    S = mf.get_ovlp()
    act = np.zeros((nao, ncas))
    for k, i in enumerate(sel):
        act[i, k] = 1.0
    # out of the core span, then Loewdin in the window's complement
    P_core = c_core @ c_core.T @ S   # the metric projector
    w = act - P_core @ act
    g = w.T @ S @ w
    gv, U = np.linalg.eigh(g)
    c_win = w @ U @ np.diag(gv ** -0.5) @ U.T
    C = np.hstack([c_core, c_win])
    P_all = C @ C.T @ S   # (C^T S C = I)
    rest = mf.mo_coeff - P_all @ mf.mo_coeff
    g2 = rest.T @ S @ rest
    gv2, U2 = np.linalg.eigh(g2)
    keep = gv2 > 1e-7 * gv2.max()
    c_virt = rest @ U2[:, keep] @ np.diag(gv2[keep] ** -0.5)
    if c_virt.shape[1] != nao - C.shape[1]:
        sys.exit(f"the complement came out {c_virt.shape[1]} wide, "
                 f"expected {nao - C.shape[1]}")
    mo = np.hstack([C, c_virt])
    # the production builder's product: the elementwise certification
    from atomterm.windows import build_full_mo, occupied_core_idx
    mo_prod = build_full_mo(mf, occupied_core_idx(mf, NCORE), sel,
                            c_core=c_core)
    # the certification: the core+window block elementwise (same
    # input, same Loewdin gauge); the virtual complement by SPAN
    # (any S-orthonormal basis of it is valid for the CAS optimizer)
    nc = C.shape[1]
    dev_blk = np.max(np.abs(mo[:, :nc] - mo_prod[:, :nc]))
    proj = mo_prod[:, nc:]
    resid = mo[:, nc:] - proj @ (proj.T @ S @ mo[:, nc:])
    dev_span = np.max(np.abs(resid))
    # the span band cites the observed: 1.9e-08 on this box (the
    # diffuse tail's overlap conditioning); the block equality is
    # the certificate
    if dev_blk > 1e-9 or dev_span > 1e-6:
        sys.exit(f"the visible window construction disagrees with "
                 f"the production builder (block dev {dev_blk:.2e}, "
                 f"complement span dev {dev_span:.2e})")
    print(f"  window: {ncas} AO columns {SHELLS}, visible "
          f"construction == production build_full_mo (block dev "
          f"{dev_blk:.1e}, complement span dev {dev_span:.1e})",
          flush=True)
    return mo, mo[:, NCORE:NCORE + ncas]


def free_leg(mol, mf, mo_full, v, ground_only):
    """The visible free route: SA -> readout -> PT2 -> meters."""
    t0 = time.time()
    mc = mcscf.CASSCF(mf, 12, SECTOR)
    mc.state_average_([1.0 / 17] * 17)
    mc.kernel(mo_full)
    print(f"  [sa] ensemble {mc.e_tot:.8f}, converged "
          f"{mc.converged} ({time.time() - t0:.1f} s)", flush=True)
    if not mc.converged:
        v.note("the SA kernel converged", "False -- the PySCF kernel "
               "hit its cycle cap; the readout still prices the "
               "frozen window")
    cas = mcscf.CASCI(mf, 12, SECTOR)
    cas.mo_coeff = mc.mo_coeff
    cas.canonicalization = False  # the window exactly as converged
    cas.fcisolver.nroots = 17
    # seed the readout FROM the SA kernel's own vectors (the two-pass
    # pattern's faithful form -- production's readout_mc(ci0=...);
    # unseeded, the davidson tails mix the quasi-degenerate block at
    # the 1e-5 level and the meters carry it)
    from pyscf.fci import cistring
    nsa = cistring.num_strings(cas.ncas, SECTOR[0])
    nsb = cistring.num_strings(cas.ncas, SECTOR[1])
    guess = [np.asarray(c).ravel().reshape(nsa, nsb) for c in mc.ci]
    cas.kernel(ci0=guess)
    expect = meters(mol, cas.mo_coeff[:, NCORE:NCORE + 12], *SECTOR)
    a0 = anchor(ADDR_412)["roots"][0]
    # per-root totals in the analyzation seam's convention
    # (e_cas[k] + the common core part -- e_tot's first entry minus
    # e_cas's first entry)
    e_tot_0 = float(np.atleast_1d(np.asarray(cas.e_tot)).ravel()[0])
    e_cas = np.atleast_1d(np.asarray(cas.e_cas)).ravel()
    core_e = e_tot_0 - float(e_cas[0])
    rows = []
    for k in (range(1) if ground_only else range(17)):
        pt = mrpt.NEVPT(cas, root=k)
        pt.ci = list(pt.ci)   # detach: the kernel writes the shared list
        e2 = float(pt.kernel())
        rows.append({"root": k, "e_ref": float(e_cas[k] + core_e),
                     "e2": e2,
                     "l2": expect("l2", cas.ci[k]),
                     "s2": expect("s2", cas.ci[k])})
        print(f"  root {k}: e_ref {rows[-1]['e_ref']:.8f}, E(2) "
              f"{e2 * 1e3:+.4f} mEh, <L2> {rows[-1]['l2']:.6f}, "
              f"<S2> {rows[-1]['s2']:.6f}", flush=True)
    r0 = rows[0]
    d_meh = (r0["e2"] - a0["e_corr"]) * 1e3
    v.check("(a) E(2) root 0 ~ the committed payload row",
            abs(d_meh) <= BAND_E2_MEH,
            f"fresh {r0['e2'] * 1e3:+.4f} mEh vs anchor "
            f"{a0['e_corr'] * 1e3:+.4f} (delta {d_meh:+.4f}, band "
            f"+-{BAND_E2_MEH})")
    dev_m = max(abs(r["l2"] - L_L2) for r in rows)
    dev_s = max(abs(r["s2"] - S_S2) for r in rows)
    v.check("(b) <L2> = 72 / <S2> = 6 on every root",
            max(dev_m, dev_s) <= BAND_L2,
            f"{len(rows)} roots, max dev {max(dev_m, dev_s):.2e} "
            f"(band {BAND_L2:g})")
    drop_row({
        "leg": "4_12", "route": "free", "cas": "CAS(4,12)",
        "e_ref": r0["e_ref"], "e2_meh": r0["e2"] * 1e3,
        "e_total": r0["e_ref"] + r0["e2"],
        "anchor_e2_meh": a0["e_corr"] * 1e3, "delta_meh": d_meh,
        "l2": f"{min(r['l2'] for r in rows):.6f}.."
              f"{max(r['l2'] for r in rows):.6f}",
        "s2": f"{max(abs(r['s2'] - S_S2) for r in rows):.1e} dev",
        "date": date.today().isoformat(),
    })
    return mc


def pinned_leg(mol, mf, v, ground_only):
    """The D2h-pinned route -- the PRODUCTION evaluator, under test:
    the lane substrate, the lane core, run_irrep_sa's two-pass
    pattern; PT2 through the per-block pinned solvers."""
    from atomterm.cases import campaign_spec, lane_substrate
    from atomterm.casscf import run_irrep_sa
    from atomterm.cells import CellKey, lane_core, window_columns
    from atomterm.windows import build_full_mo

    spec = campaign_spec("U14e")
    lmol, lmf = lane_substrate(spec)
    key = CellKey.from_address(
        "u14e_5f6d@6s6p7s_4-0_sa-casscf_ao_-_pur_avg")
    _core, nf, c_core, _rep = lane_core(key, lmol, lmf,
                                        spec.shell_bands)
    act, ncas, _w = window_columns(key, lmol, lmf, spec.shell_bands,
                                   spec.shell_radials)
    mo = build_full_mo(lmf, list(range(nf)), act, c_core=c_core)
    out = run_irrep_sa(lmol, lmf, mo[:, nf:nf + ncas], nf, *SECTOR,
                       8, "u", c_core=c_core)
    print(f"  [pinned] ensemble {out['e_tot']:.8f}, converged "
          f"{out['converged']}", flush=True)
    c_read = np.asarray(out["mo"])[:, nf:nf + ncas]
    expect = meters(lmol, c_read, *SECTOR)
    e2s, l2s = [], []
    from atomterm.core import closed_shell_core_energy
    from atomterm.irreps import branch_d2h, orb_irreps
    from atomterm.casscf import irrep_solvers
    from atomterm.perturb import nevpt2_rows_on
    from pyscf import mcscf
    orbsyms = orb_irreps(lmol, c_read)
    pinned = {irr: s for irr, s, _ in
              irrep_solvers(lmol, orbsyms, 8, "u", 1e-12)}
    core_part = closed_shell_core_energy(lmf, c_core)
    for irrep, mult in branch_d2h(8, "u").items():
        e_arr, ci = out["roots"][irrep]
        e_arr = np.atleast_1d(np.asarray(e_arr, dtype=float))
        mc = mcscf.CASCI(lmf, ncas, SECTOR)
        mc.mo_coeff = np.asarray(out["mo"], dtype=float)
        mc.fcisolver = pinned[irrep]
        mc.orbsym = orbsyms
        mc.fcisolver.orbsym = np.asarray(orbsyms)
        ci_arr = np.asarray(ci)
        mc.ci = ([ci_arr] if ci_arr.ndim == 2
                 else [np.asarray(x) for x in ci_arr])
        mc.e_cas = e_arr[:mult]
        mc.e_tot = float(mc.e_cas[0] + core_part)
        for row in nevpt2_rows_on(mc, roots=[0] if ground_only
                                  else None):
            e2s.append(row.e_corr)
            ci_vec = np.asarray(mc.ci[row.root]).ravel()
            l2s.append(expect("l2", ci_vec))
    e2s = np.array(e2s)
    print(f"  term E(2): mean {e2s.mean() * 1e3:+.4f} mEh, ptp "
          f"{np.ptp(e2s) * 1e3:.2e} mEh over {len(e2s)} components; "
          f"<L2> in [{min(l2s):.6f}, {max(l2s):.6f}]", flush=True)
    free_ens = v.free_ensemble
    delta = out["e_tot"] - free_ens
    v.check("(d) the pinned route ~ the visible free route "
            "(ensemble energies)", abs(delta) <= BAND_ROUTE,
            f"pinned {out['e_tot']:.8f} vs free {free_ens:.8f} "
            f"(delta {delta:+.2e} Eh, band {BAND_ROUTE:g})")
    drop_row({
        "leg": "4_12", "route": "pinned", "cas": "CAS(4,12)",
        "e_ref": float(out["e_tot"]), "e2_meh": float(e2s.mean() * 1e3),
        "e_total": float(out["e_tot"] + e2s.mean()),
        "anchor_e2_meh": None, "delta_meh": float(delta * 1e3),
        "l2": f"{min(l2s):.6f}..{max(l2s):.6f}",
        "s2": "72-block exact", "date": date.today().isoformat(),
    })


def mrpt_rows(mc, ground_only):
    from atomterm.perturb import nevpt2_rows_on
    return nevpt2_rows_on(mc, roots=[0] if ground_only else None)


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    ground_only = "--ground-only" in argv
    pinned = "--pinned" in argv
    v = Verdict()
    mol, mf = load_reference()
    # the substrate (marked: substrate chemistry, not under test)
    from atomterm.core import purified_core_orbitals
    c_core = purified_core_orbitals(mol, mf, NCORE)[0]
    mo_full, _c_act = visible_window(mol, mf, c_core)
    mc = free_leg(mol, mf, mo_full, v, ground_only)
    v.free_ensemble = float(mc.e_tot)
    if pinned:
        pinned_leg(mol, mf, v, ground_only)
    ok = v.report()
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
