"""The CAS(4,12) V&V leg -- the 5L leg, four routes, PySCF-direct.

THE SETUP (all visible in ``main``): the trusted ROHF rung loaded;
the core = the rotation-closed span over the occupied prefix
(atomterm.core -- substrate chemistry, marked); the window = the AO
columns of the 10f (7) + 13d (5) shells, written out and certified
elementwise against the production builder (block dev 2.2e-16).
THE FIVE ROUTES:
  free CASCI     the frozen input window, no optimizer -- the
                 pricing class (what the perturbation sees as built)
  free CASSCF    the visible window, PySCF's plain single-root kernel
  free SA        state_average_ over the ground term's 17 M_L
                 components (L = 8, u), then the two-pass readout
                 (seeded from the SA kernel's own vectors)
  pinned CI      the production lane construction (marked: under
                 test) -- per-block wfnsym solvers, no optimizer
  pinned SA      run_irrep_sa's two-pass ensemble, per-block readout
PT2: the engine's own mrpt.NEVPT per root; the meters: atomterm's
exact angular operators.

THE ASSERTIONS:
  (a0) the frozen-window CASCI ground reads the certified 5I
       (<L2> = 42, definite L) -- the SA promotes it to the 5L;
  (a) fresh SA E(2) root 0 ~ the committed payload row, band +-2 mEh;
  (b) <L2> = 72 / <S2> = 6 on every SA root, band 1e-5;
  (d) the pinned SA ensemble ~ the free SA ensemble, band 1e-6 Eh.
The free CASSCF row carries NO assertion -- the single-root draw is
the escape class by physics (the meters disclose it).

Run from the repo root:
    OMP_NUM_THREADS=16 python scripts/u14e/u14e_vv_4_12.py \
        [--ground-only]  |  tee logs/u14e_vv_4_12.log
"""

import sys
import time

import numpy as np
from pyscf import mcscf

from _vv_common import (BAND_E2_MEH, BAND_L2, BAND_ROUTE, L_L2, S_S2,
                        Verdict, anchor, conv_of, core_part_of, drop_row,
                        fmt_range, free_kernel, lane_route,
                        load_reference, make_blocks, meters, meters_all,
                        nevpt_rows, pinned_ci_run, pinned_sa_run, price,
                        std_row)

ADDR_412 = "u14e_10f13d_4-0_sa-casscf_ao_-_pur_avg"
ADDR_412_LANE = "u14e_5f6d@6s6p7s_4-0_sa-casscf_ao_-_pur_avg"
ADDR_412_CASSCF = "u14e_10f13d_4-0_casscf_ao_-_pur_avg"
SECTOR = (4, 0)
NCORE = 5
SHELLS = (("10f", 7), ("13d", 5))
L, PARITY = 8, "u"
L2_5I = 42.0  # the certified INPUT-window ground: L = 6 -> L(L+1)
ROOTS = 17


def visible_window(mol, mf, c_core):
    """The window, written out: select the shells' AO columns, project
    the core span out (the S metric), Loewdin-orthonormalize, fill the
    complement from the canonical virtuals.  Certified elementwise
    (core+window) against the production builder; the complement by
    span."""
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
    from atomterm.windows import build_full_mo, occupied_core_idx
    mo_prod = build_full_mo(mf, occupied_core_idx(mf, NCORE), sel,
                            c_core=c_core)
    nc = C.shape[1]
    dev_blk = np.max(np.abs(mo[:, :nc] - mo_prod[:, :nc]))
    proj = mo_prod[:, nc:]
    resid = mo[:, nc:] - proj @ (proj.T @ S @ mo[:, nc:])
    dev_span = np.max(np.abs(resid))
    # the span band cites the observed: 1.9e-08 (the diffuse tail's
    # overlap conditioning); the block equality is the certificate
    if dev_blk > 1e-9 or dev_span > 1e-6:
        sys.exit(f"the visible window construction disagrees with "
                 f"the production builder (block dev {dev_blk:.2e}, "
                 f"complement span dev {dev_span:.2e})")
    print(f"  window: {ncas} AO columns {SHELLS}, visible "
          f"construction == production build_full_mo (block dev "
          f"{dev_blk:.1e}, complement span dev {dev_span:.1e})",
          flush=True)
    return mo


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    ground_only = "--ground-only" in argv
    pt2_roots = [0] if ground_only else None
    v = Verdict()
    mol, mf = load_reference()
    from atomterm.core import purified_core_orbitals
    c_core = purified_core_orbitals(mol, mf, NCORE)[0]
    mo_full = visible_window(mol, mf, c_core)
    expect = meters(mol, mo_full[:, NCORE:NCORE + 12], *SECTOR)
    weights = [1.0 / 17] * 17

    # ---- free CASCI (the frozen input window; the pricing class) -----
    # no optimizer, deterministic: the window prices the term AS BUILT
    # (the campaign's "fuel" reading -- what the perturbation sees
    # before any relaxation)
    t0 = time.time()
    cas_f = mcscf.CASCI(mf, 12, SECTOR)
    cas_f.canonicalization = False
    cas_f.kernel(mo_full)
    r_f = price(cas_f, expect, roots=[0])[0]
    print(f"  [casci] e_tot {cas_f.e_tot:.8f}, <L2> {r_f['l2']:.6f} "
          f"({time.time() - t0:.1f} s)", flush=True)
    v.check("(a0) the frozen-window CASCI ground reads the certified "
            "5I (<L2> = 42, definite L -- the input window is "
            "rotation-closed, the (4,0) sector complete; the SA then "
            "promotes it to the 5L -- the certified record's 'the ~1 "
            "Eh gain IS the 5L relaxation')",
            abs(r_f["l2"] - L2_5I) <= BAND_L2,
            f"<L2> {r_f['l2']:.6f} vs 42 (band {BAND_L2:g})")
    drop_row(std_row("casci", "free@pur", "(4,12)", float(cas_f.e_tot),
                     r_f["e2"] * 1e3,
                     float(cas_f.e_tot) + r_f["e2"],
                     f"{r_f['l2']:.4f}", f"{r_f['s2']:.4f}",
                     "5I input ground",
                     conv=conv_of(cas_f.fcisolver.conv_tol, None)))

    # ---- free CASSCF (single root; the escape class, no assertion) --
    t0 = time.time()
    mc_c, cas_c = free_kernel(mf, 12, SECTOR, mo_full, track_iters=True)
    r_c = price(cas_c, expect, roots=[0])[0]
    a_c = anchor(ADDR_412_CASSCF)["roots"][0]
    print(f"  [casscf] e_tot {mc_c.e_tot:.8f}, converged "
          f"{mc_c.converged}, <L2> {r_c['l2']:.6f} "
          f"({time.time() - t0:.1f} s)", flush=True)
    v.note("free CASSCF (the single-root draw class)",
           f"<L2> {r_c['l2']:.6f}, E(2) {r_c['e2'] * 1e3:+.4f} mEh "
           f"(the committed draw priced {a_c['e_corr'] * 1e3:+.4f}; "
           f"the single-root endpoints differ by draw -- the meters "
           f"are the disclosure)")
    drop_row(std_row("casscf", "free@pur", "(4,12)", float(mc_c.e_tot),
                     r_c["e2"] * 1e3,
                     float(mc_c.e_tot) + r_c["e2"],
                     f"{r_c['l2']:.4f}", f"{r_c['s2']:.4f}",
                     "escape draw",
                     conv=conv_of(mc_c.final_de, mc_c.n_macro,
                                  bool(mc_c.converged))))

    # ---- free SA ------------------------------------------------------
    t0 = time.time()
    mc_sa, cas_sa = free_kernel(mf, 12, SECTOR, mo_full,
                                weights=weights, nroots=ROOTS,
                                track_iters=True)
    core_e, e_cas = core_part_of(cas_sa)
    expect_sa = meters(mol, cas_sa.mo_coeff[:, NCORE:NCORE + 12],
                       *SECTOR)
    a0 = anchor(ADDR_412)["roots"][0]
    rows = price(cas_sa, expect_sa,
                 roots=range(1) if ground_only else range(ROOTS))
    for r in rows:
        r["e_ref"] = float(e_cas[r["root"]] + core_e)
        print(f"  root {r['root']}: e_ref {r['e_ref']:.8f}, E(2) "
              f"{r['e2'] * 1e3:+.4f} mEh, <L2> {r['l2']:.6f}, "
              f"<S2> {r['s2']:.6f}", flush=True)
    r0 = rows[0]
    d_meh = (r0["e2"] - a0["e_corr"]) * 1e3
    v.check("(a) SA E(2) root 0 ~ the committed payload row",
            abs(d_meh) <= BAND_E2_MEH,
            f"fresh {r0['e2'] * 1e3:+.4f} mEh vs anchor "
            f"{a0['e_corr'] * 1e3:+.4f} (delta {d_meh:+.4f}, band "
            f"+-{BAND_E2_MEH})")
    dev_m = max(abs(r["l2"] - L_L2) for r in rows)
    dev_s = max(abs(r["s2"] - S_S2) for r in rows)
    v.check("(b) <L2> = 72 / <S2> = 6 on every SA root",
            max(dev_m, dev_s) <= BAND_L2,
            f"{len(rows)} roots, max dev {max(dev_m, dev_s):.2e} "
            f"(band {BAND_L2:g})")
    drop_row(std_row("sa-casscf", "free@pur", "(4,12)",
                     float(e_cas[0] + core_e), r0["e2"] * 1e3,
                     float(r0["e_ref"] + r0["e2"]),
                     fmt_range([r['l2'] for r in rows]),
                     f"{max(r['s2'] for r in rows):.4f}", "5L",
                     conv=conv_of(mc_sa.final_de, mc_sa.n_macro,
                                  bool(mc_sa.converged))))
    print(f"  [sa] ensemble {mc_sa.e_tot:.8f}, converged "
          f"{mc_sa.converged} ({time.time() - t0:.1f} s)", flush=True)

    # ---- the pinned substrate (built once, two legs) ------------------
    from atomterm.cases import campaign_spec
    lmol, lmf, c_core_l, nf, mo_l, ncas_l, _act = lane_route(
        campaign_spec("U14e"), ADDR_412_LANE, SECTOR)

    # ---- pinned CI (per-block solvers, no optimizer) ------------------
    t0 = time.time()
    roots_ci, _ = pinned_ci_run(lmol, lmf, mo_l, nf, ncas_l, SECTOR,
                                L, PARITY, c_core_l)
    blocks, expect_p, _ = make_blocks(lmol, lmf, mo_l, nf, ncas_l,
                                      SECTOR, L, PARITY, c_core_l,
                                      roots_ci)
    l2s, s2s = meters_all(blocks, expect_p)
    e2s = [row.e_corr for _, mc in blocks
           for row in nevpt_rows(mc, pt2_roots)]
    term_means = []
    for _, mc in blocks:
        ec = np.atleast_1d(mc.e_cas).ravel()
        cb = float(mc.e_tot) - float(ec[0])
        term_means.append(float(np.mean(ec + cb)))
    e_var_ci = float(np.mean(term_means))
    print(f"  [pinned casci] term e_var {e_var_ci:.8f}, E(2) mean "
          f"{np.mean(e2s) * 1e3:+.4f} mEh, ptp "
          f"{np.ptp(e2s) * 1e3:.2e} mEh, <L2> "
          f"[{min(l2s):.6f}, {max(l2s):.6f}] "
          f"({time.time() - t0:.1f} s)", flush=True)
    v.note("pinned casci (the pinned analog of the unoptimized row)",
           f"<L2> in [{min(l2s):.4f}, {max(l2s):.4f}], <S2> max "
           f"{max(s2s):.4f}; E(2) mean {np.mean(e2s) * 1e3:+.4f} mEh")
    # the lane-core route is outside the table's keep-set (the
    # pinned@pur rows are the table's pinned route): disclosed only
    v.note("the lane-pinned casci (print-only)",
           f"E_method {e_var_ci:.8f}, E(2) mean "
           f"{np.mean(e2s) * 1e3:+.4f} mEh")

    # ---- pinned SA ----------------------------------------------------
    t0 = time.time()
    out_sa, _ = pinned_sa_run(lmol, lmf, mo_l, nf, ncas_l, SECTOR,
                              L, PARITY, c_core_l)
    blocks_s, expect_s, _ = make_blocks(lmol, lmf,
                                        np.asarray(out_sa["mo"]),
                                        nf, ncas_l, SECTOR, L, PARITY,
                                        c_core_l, out_sa["roots"])
    l2s_s, s2s_s = meters_all(blocks_s, expect_s)
    e2s_s = [row.e_corr for _, mc in blocks_s
             for row in nevpt_rows(mc, pt2_roots)]
    print(f"  [pinned sa] ensemble {out_sa['e_tot']:.8f}, term E(2) "
          f"mean {np.mean(e2s_s) * 1e3:+.4f} mEh, ptp "
          f"{np.ptp(e2s_s) * 1e3:.2e} mEh, <L2> "
          f"[{min(l2s_s):.6f}, {max(l2s_s):.6f}] "
          f"({time.time() - t0:.1f} s)", flush=True)
    delta = float(out_sa["e_tot"]) - float(mc_sa.e_tot)
    v.check("(d) the pinned SA ensemble ~ the free SA ensemble",
            abs(delta) <= BAND_ROUTE,
            f"pinned {out_sa['e_tot']:.8f} vs free {mc_sa.e_tot:.8f} "
            f"(delta {delta:+.2e} Eh, band {BAND_ROUTE:g})")
    v.note("the lane-pinned SA (print-only)",
           f"ensemble {out_sa['e_tot']:.8f}, E(2) mean "
           f"{np.mean(e2s_s) * 1e3:+.4f} mEh")

    ok = v.report()
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
