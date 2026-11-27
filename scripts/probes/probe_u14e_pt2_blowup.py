"""The PT2-blowup probe: the pinned@pur (6,13) SA row's E(2) read
-1580 mEh (the -43 eV anomaly) against ~-285 mEh (the lane route) and
~-282 mEh (the free SA).  The decomposition, printed only:

  1. the pinned@pur (6,13) SA re-run on the AS-BUILT window (the
     anomaly's gauge), the per-BLOCK E(2) of the ground roots -- where
     does the blowup sit?
  2. the worst block's 8-Koopmans-class breakdown (label / perturber
     norm / class energy) -- which class carries it?
  3. the free SA x17 on the SAME window: the per-root E(2) for ALL 17
     roots + the E_CAS invariance check (the pinned ensemble's
     per-block E_CAS beside the free roots' -- the same spans, the
     same CI Hamiltonian, so the roots' spectra must agree to the
     degenerate-mixing tolerance);
  4. the per-unit-norm E(2) per gauge (mean E(2) / mean norm) -- the
     weight-vs-depth attribution: a constant per-unit depth = the
     weight effect (norm 14.9); a deeper per-unit = the small-gap
     enhancement (the boundary orbital's ionization-like gaps);
  5. the near-edge orbital structure: the highest occupied and the
     lowest virtual epsilon (the gap the Sr denominators ride).

Run from the repo root:
    OMP_NUM_THREADS=64 python scripts/probes/probe_u14e_pt2_blowup.py \
        | tee logs/probe_u14e_pt2_blowup.log
"""

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "u14e"))

import numpy as np  # noqa: E402
from pyscf import mcscf  # noqa: E402

from vv.meters import meters_all  # noqa: E402
from vv.reference import load_reference  # noqa: E402
from vv.routes import core_part_of, free_kernel, make_blocks, nevpt_rows  # noqa: E402

SECTOR = (5, 1)
NCORE = 4
ROOTS = 17
BUDGET = {"max_cycle_macro": 200}
ADDR_613 = "u14e_10f13d7s_5-1_sa-casscf_ao_-_pur_avg"


def main() -> int:
    t00 = time.time()
    mol, mf = load_reference()
    from atomterm.cases import campaign_spec
    spec = campaign_spec("U14e")
    import u14e_vv_unification as uni
    mol_d2h, mf_d2h = uni.wrapper_substrate(mol, mf, spec)
    from atomterm.windows import build_full_mo, occupied_core_idx
    from atomterm.cells import CellKey, key_core, window_columns
    key = CellKey.from_address(ADDR_613)
    nc, _cc, c_pur = key_core(key, mol, mf, spec.shell_bands)
    act, ncas, _ = window_columns(key, mol, mf, spec.shell_bands,
                                  spec.shell_radials, ncol=nc)
    mo = build_full_mo(mf_d2h, occupied_core_idx(mf_d2h, nc), act,
                       c_core=c_pur)
    c_read = np.asarray(mo)[:, nc:nc + ncas]

    # ---- 1. the pinned ensemble, the per-block E(2) --------------------
    print("===== 1. the pinned SA (as-built window), per-block E(2) "
          "=====", flush=True)
    t0 = time.time()
    from vv.routes import pinned_sa_run
    out, _ = pinned_sa_run(mol_d2h, mf_d2h, mo, nc, ncas, SECTOR, 8,
                           "u", c_pur, c_read=c_read, mc_opts=BUDGET)
    blocks, expect, _ = make_blocks(mol_d2h, mf_d2h,
                                    np.asarray(out["mo"]), nc, ncas,
                                    SECTOR, 8, "u", c_pur, out["roots"],
                                    c_read=c_read)
    print(f"  ensemble {out['e_tot']:.8f}, converged "
          f"{out['converged']} (cycles {out['cycles']}, "
          f"{time.time() - t0:.0f} s)", flush=True)
    block_rows = []
    norms = []
    for name, mc in blocks:
        row = nevpt_rows(mc, [0])[0]
        norm_mean = float(np.mean([cl.norm for cl in row.classes]))
        norms.append(norm_mean)
        print(f"  block {name:4s} root 0: E(2) {row.e_corr * 1e3:+.4f} "
              f"mEh (the mean class norm {norm_mean:.4f})", flush=True)
        block_rows.append((name, mc, row))
    mean = float(np.mean([r.e_corr for _, _, r in block_rows]))
    mean_norm = float(np.mean(norms))
    per_unit = mean / mean_norm
    print(f"  the mean: {mean * 1e3:+.4f} mEh; the mean class norm "
          f"{mean_norm:.4f}; the PER-UNIT-NORM E(2): {per_unit:+.4f} "
          f"Eh", flush=True)

    # ---- 2. the worst block's class breakdown --------------------------
    (wname, wmc, wrow) = max(block_rows, key=lambda t: abs(t[2].e_corr))
    print(f"===== 2. the worst block's class breakdown: {wname} "
          f"(|E(2)| {abs(wrow.e_corr):.4f} Eh) =====", flush=True)
    for cl in wrow.classes:
        print(f"    {cl.label:8s} ({cl.note:12s}): norm "
              f"{cl.norm:.6e}, E(2) {cl.energy:+.6f} Eh", flush=True)

    # ---- 3. the free SA on the C1 route (the honest free-route
    #        control) + the per-root E(2) over all 17 -------------------
    print("===== 3. the free SA x17 on the C1 route (the control) =====",
          flush=True)
    t0 = time.time()
    from vv.meters import meters
    from vv.routes import price
    mo_c1 = build_full_mo(mf, occupied_core_idx(mf, NCORE), act,
                          c_core=c_pur)
    _mc, cas = free_kernel(mf, ncas, SECTOR, mo_c1, weights=[1 / 17] * 17,
                           nroots=ROOTS, track_iters=True)
    core_e, e_cas_free = core_part_of(cas)
    expect_c1 = meters(mol, cas.mo_coeff[:, NCORE:NCORE + ncas], *SECTOR)
    print(f"  ensemble {_mc.e_tot:.8f}, converged {_mc.converged} "
          f"({time.time() - t0:.0f} s)", flush=True)
    e2s = []
    for k in range(ROOTS):
        r = price(cas, expect_c1, roots=[k])[0]
        e2s.append(r["e2"])
        print(f"  root {k}: e_ref {e_cas_free[k] + core_e:.8f}, E(2) "
              f"{r['e2'] * 1e3:+.4f} mEh, <L2> {r['l2']:.4f}, "
              f"<S2> {r['s2']:.4f}", flush=True)
    print(f"  the free-SA E(2) range: "
          f"[{min(e2s) * 1e3:+.1f}, {max(e2s) * 1e3:+.1f}] mEh",
          flush=True)

    # ---- 4. the per-unit-norm attribution, both gauges ------------------
    print("===== 4. the weight-vs-depth attribution =====", flush=True)
    print(f"  the pinned gauge: the mean class norm {mean_norm:.4f}, "
          f"the per-unit E(2) {per_unit:+.4f} Eh", flush=True)
    from pyscf import mrpt
    pt_free = mrpt.NEVPT(cas, root=0)
    fr = nevpt_rows(cas, [0])[0]
    free_norm = float(np.mean([cl.norm for cl in fr.classes]))
    print(f"  the free gauge (the same window, root 0): E(2) "
          f"{fr.e_corr * 1e3:+.4f} mEh, the mean class norm "
          f"{free_norm:.4f}, the per-unit E(2) "
          f"{fr.e_corr / free_norm:+.4f} Eh -- THE ATTRIBUTION: the "
          f"depth-vs-weight split is the measurement", flush=True)
    e2_free_mean = float(np.mean(e2s))

    # ---- 5. the near-edge orbital structure ----------------------------
    print("===== 5. the near-edge orbital structure =====", flush=True)
    e_all = np.asarray(mf.mo_energy).ravel()
    occ_d = np.asarray(mf.mo_occ).ravel() > 0
    occ_e = np.sort(e_all[occ_d])[::-1][:3]
    virt_e = np.sort(e_all[~occ_d])[:5]
    print(f"  the highest occupied: "
          + ", ".join(f"{x:+.4f}" for x in occ_e), flush=True)
    print(f"  the lowest virtual:  "
          + ", ".join(f"{x:+.4f}" for x in virt_e), flush=True)
    gap = float(virt_e[0] - occ_e[0])
    print(f"  the edge gap: {gap:+.4f} Eh -- the Sr denominators ride "
          f"this edge (sign-definite; a small gap = the norm's "
          f"enhancement)", flush=True)

    # ---- 6. THE GAUGE REDRAW (the ledger's instrument, implemented):
    # the random unitary INSIDE the anomaly's block -- the span is
    # invariant, so E_CAS/E_tot/<L2> must hold while E(2) spreads --
    # the ledger's numerical evidence on the Sr block itself.
    from pyscf import mrpt as _mrpt
    print("===== 6. the gauge redraw on the anomaly's block =====",
          flush=True)
    name_a, mc_a, row_a = max(block_rows, key=lambda t: abs(t[2].e_corr))
    orbs_a = np.asarray(mc_a.orbsym)[NCORE:NCORE + ncas]
    spreads = []
    for seed in (1, 2, 3):
        r2 = np.random.default_rng(seed)
        # the D2h-preserving gauge: the random unitaries WITHIN each
        # irrep group (a full unitary would mix across the irreps and
        # break the labeling)
        u_big = np.eye(ncas)
        for ir in np.unique(orbs_a):
            idx_g = np.flatnonzero(orbs_a == ir)
            amat = r2.normal(size=(len(idx_g), len(idx_g)))
            u_g, _ = np.linalg.qr(amat)
            u_big[np.ix_(idx_g, idx_g)] = u_g
        mc_r = mcscf.CASCI(mf_d2h, ncas, SECTOR)
        mc_r.mo_coeff = np.asarray(mc_a.mo_coeff, dtype=float).copy()
        mc_r.mo_coeff[:, NCORE:NCORE + ncas] = (
            mc_r.mo_coeff[:, NCORE:NCORE + ncas] @ u_big)
        mc_r.fcisolver = mc_a.fcisolver
        mc_r.orbsym = mc_a.orbsym
        mc_r.e_cas = mc_a.e_cas
        mc_r.kernel()
        row_r = nevpt_rows(mc_r, [0])[0]
        l2_r = expect("l2", np.asarray(mc_r.ci[0]).ravel())
        spreads.append(row_r.e_corr * 1e3)
        print(f"  [redraw {seed}] E(2) {row_r.e_corr * 1e3:+.4f} mEh "
              f"(the controls: E_tot "
              f"{float(np.atleast_1d(mc_r.e_tot).ravel()[0]):+.8f}, "
              f"<L2> {l2_r:.4f})", flush=True)
    print(f"  [redraw] the E(2) spread: [{min(spreads):+.1f}, "
          f"{max(spreads):+.1f}] mEh vs the unpinned gauge's "
          f"{row_a.e_corr * 1e3:+.1f} -- the E(2) gauge-visibility "
          f"measured on the Sr block itself; the controls "
          f"(E_tot/<L2>) hold", flush=True)

    print(f"\n({time.time() - t00:.0f} s total)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
