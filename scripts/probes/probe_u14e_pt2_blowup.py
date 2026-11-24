"""The PT2-blowup probe: the pinned@pur (6,13) SA row's E(2) read
-1580 mEh (the -43 eV anomaly) against ~-285 mEh (the lane route) and
~-282 mEh (the free SA).  The decomposition, printed only:

  1. the pinned@pur (6,13) SA re-run on the AS-BUILT window (the
     anomaly's gauge), the per-BLOCK E(2) of the ground roots -- where
     does the blowup sit?
  2. the worst block's 8-Koopmans-class breakdown (label / perturber
     norm / class energy) -- which class diverges?
  3. the free SA x17 on the SAME window, the per-root E(2) -- the
     normal contrast (isolating the pinned-block plumbing vs the
     gauge).

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

from _vv_common import (core_part_of, free_kernel, load_reference,  # noqa: E402
                        make_blocks, nevpt_rows)

SECTOR_613 = (5, 1)
NCORE_613 = 4
SECTOR_412 = (4, 0)
NCORE_412 = 5
ROOTS = 17
BUDGET = {"max_cycle_macro": 200}
ADDR_613 = "u14e_10f13d7s_5-1_sa-casscf_ao_-_pur_avg"
ADDR_412 = "u14e_10f13d_4-0_casci_ao_-_pur_avg"


def _decompose(tag, mol_d2h, mf_d2h, mo, nc, ncas, c_pur, c_read,
               sector):
    """The pinned ensemble on one window: the SA, the per-block E(2)
    of the ground roots, the worst block's class breakdown."""
    t0 = time.time()
    from _vv_common import pinned_sa_run
    out, _ = pinned_sa_run(mol_d2h, mf_d2h, mo, nc, ncas, sector,
                           8, "u", c_pur, c_read=c_read, mc_opts=BUDGET)
    blocks, expect, _ = make_blocks(mol_d2h, mf_d2h,
                                    np.asarray(out["mo"]), nc, ncas,
                                    sector, 8, "u", c_pur, out["roots"],
                                    c_read=c_read)
    print(f"  [{tag}] ensemble {out['e_tot']:.8f}, converged "
          f"{out['converged']} (cycles {out['cycles']}, "
          f"{time.time() - t0:.0f} s)", flush=True)
    block_rows = []
    for name, mc in blocks:
        row = nevpt_rows(mc, [0])[0]
        print(f"  [{tag}] block {name:4s} root 0: E(2) "
              f"{row.e_corr * 1e3:+.4f} mEh", flush=True)
        block_rows.append((name, mc, row))
    mean = float(np.mean([r.e_corr for _, _, r in block_rows]))
    print(f"  [{tag}] the mean: {mean * 1e3:+.4f} mEh", flush=True)
    (wname, wmc, wrow) = max(block_rows,
                             key=lambda t: abs(t[2].e_corr))
    print(f"  [{tag}] the worst block's classes: {wname}", flush=True)
    for cl in wrow.classes:
        print(f"    {cl.label:8s} ({cl.note:12s}): norm "
              f"{cl.norm:.6e}, E(2) {cl.energy:+.6f} Eh", flush=True)


def main() -> int:
    t00 = time.time()
    mol, mf = load_reference()
    from atomterm.cases import campaign_spec
    spec = campaign_spec("U14e")
    import u14e_vv_unification as uni
    mol_d2h, mf_d2h = uni.wrapper_substrate(mol, mf, spec)
    from atomterm.windows import build_full_mo, occupied_core_idx
    from atomterm.cells import CellKey, key_core, window_columns

    # ---- case 1: the (6,13) as-built window ----------------------------
    print("===== case 1: the (6,13) as-built window =====", flush=True)
    key = CellKey.from_address(ADDR_613)
    nc, _cc, c_pur = key_core(key, mol, mf, spec.shell_bands)
    act, ncas, _ = window_columns(key, mol, mf, spec.shell_bands,
                                  spec.shell_radials, ncol=nc)
    mo = build_full_mo(mf_d2h, occupied_core_idx(mf_d2h, nc), act,
                       c_core=c_pur)
    c_read = np.asarray(mo)[:, nc:nc + ncas]
    _decompose("(6,13)", mol_d2h, mf_d2h, mo, nc, ncas, c_pur, c_read,
               SECTOR_613)

    # ---- case 2: the (4,12) raw window ---------------------------------
    print("===== case 2: the (4,12) raw window =====", flush=True)
    key2 = CellKey.from_address(ADDR_412)
    nc2, _cc2, c_pur2 = key_core(key2, mol, mf, spec.shell_bands)
    act2, ncas2, _ = window_columns(key2, mol, mf, spec.shell_bands,
                                    spec.shell_radials, ncol=nc2)
    mo2 = build_full_mo(mf_d2h, occupied_core_idx(mf_d2h, nc2), act2,
                        c_core=c_pur2)
    c_read2 = np.asarray(mo2)[:, nc2:nc2 + ncas2]
    _decompose("(4,12)", mol_d2h, mf_d2h, mo2, nc2, ncas2, c_pur2,
               c_read2, SECTOR_412)

    # ---- the contrast: the free SA on the same windows -----------------
    print("===== the contrast: the free SA per-root E(2) =====",
          flush=True)
    from _vv_common import price, meters
    for tag, mo_w, ncw, ncasw, sector in (
            ("(6,13)", mo, nc, ncas, SECTOR_613),
            ("(4,12)", mo2, nc2, ncas2, SECTOR_412)):
        _mc, cas = free_kernel(mf_d2h, ncasw, sector, mo_w,
                               weights=[1 / 17] * 17, nroots=ROOTS,
                               track_iters=True)
        expect = meters(mol_d2h, cas.mo_coeff[:, ncw:ncw + ncasw],
                        *sector)
        e2s = [r["e2"] * 1e3 for k in (0, ROOTS - 1)
               for r in [price(cas, expect, roots=[k])[0]]]
        print(f"  [{tag}] the free SA: E(2) root 0 {e2s[0]:+.1f}, "
              f"root {ROOTS - 1} {e2s[1]:+.1f} mEh", flush=True)
    print(f"({time.time() - t00:.0f} s total)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

