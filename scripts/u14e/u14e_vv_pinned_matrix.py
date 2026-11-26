"""The pinned-matrix completion leg -- the pinned@pur route's missing
rows, so the pinned routes' structure mirrors the free routes' CAS
ladder exactly (the user's directive; the rest of the matrix is dash
rows rendered by ``render_table``).

THE FIVE RUNS (the wrapper route throughout: the C1 rung's
gauge-independent AO integrals on the D2h-tagged mol, the purified
core, NO SCF):
  casci (4,12)           the RAW-band window (10f13d) through the
                         pinned route -- its block grounds are the
                         class disclosure (the free casci's raw
                         window freezes the 5I; the pinned blocks
                         read per-irrep grounds);
  casci (4,15)@r1        the prior-axis f/d radials + the promoted
                         14p band at identity;
  sa-casscf (4,15)@r1    the same window, the term-ensemble
                         optimizer;
  sa-casscf (6,13)       the AS-BUILT distilled window -- the
                         per-component <L2> walk prints the 72/90
                         counts (the user's question: how much of
                         the declared term survives the pinned
                         ensemble's optimization);
  sa-casscf (6,13)@r1    the same walk on step B's window.

The existing C1-gate rows were RE-LABELED (4,12) -> (4,12)@r1 in the
store (they ran the @r1 window; the measurements stand).

Run from the repo root:
    OMP_NUM_THREADS=64 python scripts/u14e/u14e_vv_pinned_matrix.py \
        | tee logs/u14e_vv_pinned_matrix.log
"""

import sys
import time

import numpy as np

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))

from _vv_common import (BUDGET_200, Verdict, conv_of, drop_row, fmt_range,  # noqa: E402
                        load_reference, make_blocks, meters_all,
                        nevpt_rows, pinned_ci_run, pinned_sa_run,
                        std_row)

ADDR_412_RAW = "u14e_10f13d_4-0_casci_ao_-_pur_avg"
ADDR_415_R1 = "u14e_5f6d14p@r1_4-0_casci_ao_-_pur_avg"
ADDR_613 = "u14e_10f13d7s_5-1_sa-casscf_ao_-_pur_avg"
ADDR_412_R1 = "u14e_5f6d@r1_4-0_casci_ao_-_pur_avg"
L, PARITY = 8, "u"
BUDGET = BUDGET_200  # the optimizer-budget directive (from _vv_common)


def _counts(tag, l2s):
    n72 = sum(abs(x - 72.0) <= 0.01 for x in l2s)
    n90 = sum(abs(x - 90.0) <= 0.01 for x in l2s)
    print(f"  [{tag}] the counts: 72.0 +- 0.01: {n72}; 90.0 +- 0.01: "
          f"{n90}; other: {len(l2s) - n72 - n90}", flush=True)
    return n72, n90


def _ci_section(v, tag, addr, mol, mf, spec, mol_d2h, mf_d2h, sector,
                ncore):
    """One pinned@pur CI row: the window through the wrapper, the
    block grounds + meters, PT2 on the ground roots, the row."""
    from atomterm.windows import build_full_mo, occupied_core_idx
    import u14e_vv_unification as uni
    mo, ncas, c_pur, nf = uni.window_of(addr, mol, mf, spec, mf_d2h)
    t0 = time.time()
    roots, _ = pinned_ci_run(mol_d2h, mf_d2h, mo, nf, ncas, sector,
                             L, PARITY, c_pur)
    blocks, expect, _ = make_blocks(mol_d2h, mf_d2h, mo, nf, ncas,
                                    sector, L, PARITY, c_pur, roots)
    l2s, s2s = meters_all(blocks, expect)
    print(f"  [{tag} casci] {len(l2s)} components, <L2> in "
          f"[{min(l2s):.6f}, {max(l2s):.6f}] ({time.time() - t0:.1f} s)",
          flush=True)
    _counts(tag, l2s)
    e2s = [row.e_corr for _, mc in blocks
           for row in nevpt_rows(mc, [0])]
    term_means = []
    for _, mc in blocks:
        ec = np.atleast_1d(mc.e_cas).ravel()
        cb = float(mc.e_tot) - float(ec[0])
        term_means.append(float(np.mean(ec + cb)))
    e_var_ci = float(np.mean(term_means))
    drop_row(std_row("casci", "pinned@pur", tag, e_var_ci,
                     float(np.mean(e2s) * 1e3),
                     float(e_var_ci + np.mean(e2s)),
                     fmt_range(l2s), f"{max(s2s):.4f}",
                     "block grounds",
                     conv=conv_of(1e-12, None)))


def _sa_section(v, tag, addr_or_mo, mol, mf, spec, mol_d2h, mf_d2h,
                sector, ncore, walk=False):
    """One pinned@pur SA row: the ensemble optimizer, the exact
    per-component readout, the per-component <L2> walk when ``walk``
    (the user's 72/90 counts), PT2 on the ground roots, the row."""
    from atomterm.windows import build_full_mo, occupied_core_idx
    import u14e_vv_unification as uni
    if isinstance(addr_or_mo, str):
        mo, ncas, c_pur, nf = uni.window_of(addr_or_mo, mol, mf, spec,
                                            mf_d2h)
        c_read = np.asarray(mo)[:, nf:nf + ncas]
    else:
        mo, c_pur, nf, ncas = addr_or_mo
        c_read = np.asarray(mo)[:, nf:nf + ncas]
    t0 = time.time()
    out, _ = pinned_sa_run(mol_d2h, mf_d2h, mo, nf, ncas, sector,
                           L, PARITY, c_pur, c_read=c_read,
                           mc_opts=BUDGET)
    blocks, expect, _ = make_blocks(mol_d2h, mf_d2h,
                                    np.asarray(out["mo"]), nf, ncas,
                                    sector, L, PARITY, c_pur,
                                    out["roots"], c_read=c_read)
    l2s, s2s = meters_all(blocks, expect)
    print(f"  [{tag} sa] ensemble {out['e_tot']:.8f}, converged "
          f"{out['converged']} (cycles {out['cycles']}, "
          f"{time.time() - t0:.1f} s); <L2> [{min(l2s):.6f}, "
          f"{max(l2s):.6f}]", flush=True)
    if walk:
        i = 0
        for name, mc in blocks:
            for k in range(len(np.atleast_1d(mc.e_cas))):
                ci = np.asarray(mc.ci[k]).ravel()
                print(f"    {tag} {name}[{k}]: <L2> "
                      f"{expect('l2', ci):.6f}, <S2> "
                      f"{expect('s2', ci):.6f}", flush=True)
                i += 1
    _counts(tag, l2s)
    e2s = [row.e_corr for _, mc in blocks
           for row in nevpt_rows(mc, [0])]
    drop_row(std_row("sa-casscf", "pinned@pur", tag,
                     float(out["e_tot"]),
                     float(np.mean(e2s) * 1e3),
                     float(out["e_tot"] + np.mean(e2s)),
                     fmt_range(l2s), f"{max(s2s):.4f}",
                     "draw-dependent partition",
                     conv=conv_of(out["conv_dE"], out["cycles"],
                                  bool(out["converged"]))))


def main() -> int:
    v = Verdict()
    mol, mf = load_reference()
    from atomterm.cases import campaign_spec
    spec = campaign_spec("U14e")
    import u14e_vv_unification as uni
    mol_d2h, mf_d2h = uni.wrapper_substrate(mol, mf, spec)

    # ---- casci (4,12): the RAW-band window through the pinned route
    _ci_section(v, "(4,12)", ADDR_412_RAW, mol, mf, spec, mol_d2h,
                mf_d2h, (4, 0), 5)

    # ---- casci (4,15)@r1 ----------------------------------------------
    _ci_section(v, "(4,15)@r1", ADDR_415_R1, mol, mf, spec, mol_d2h,
                mf_d2h, (4, 0), 5)

    # ---- sa-casscf (4,15)@r1 ------------------------------------------
    _sa_section(v, "(4,15)@r1", ADDR_415_R1, mol, mf, spec, mol_d2h,
                mf_d2h, (4, 0), 5)

    # ---- sa-casscf (6,13): the as-built window, the full walk ---------
    from atomterm.windows import build_full_mo, occupied_core_idx
    from atomterm.cells import CellKey, key_core, window_columns
    key = CellKey.from_address(ADDR_613)
    nc613, _cc, c_pur613 = key_core(key, mol, mf, spec.shell_bands)
    act613, ncas613, _ = window_columns(key, mol, mf, spec.shell_bands,
                                        spec.shell_radials, ncol=nc613)
    mo613 = build_full_mo(mf_d2h, occupied_core_idx(mf_d2h, nc613),
                          act613, c_core=c_pur613)
    _sa_section(v, "(6,13)", (mo613, c_pur613, nc613, ncas613), mol,
                mf, spec, mol_d2h, mf_d2h, (5, 1), nc613, walk=True)

    # ---- sa-casscf (6,13)@r1: step B's window, the full walk ----------
    act_s, ncas_s, _purity = uni.act_s_mo_of(mol, mf, mol_d2h, spec)
    mo_s = build_full_mo(mf_d2h, occupied_core_idx(mf_d2h, nc613),
                         act_s, c_core=c_pur613)
    _sa_section(v, "(6,13)@r1", (mo_s, c_pur613, nc613, ncas_s), mol,
                mf, spec, mol_d2h, mf_d2h, (5, 1), nc613, walk=True)

    ok = v.report()
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
