"""The CAS(6,13) V&V leg -- the 7s-promotion leg of the flagship
CAS-size axis, sector (5,1).  THE DISCLOSURE LEG: everything is
printed; nothing gates (the stored (6,13) SA cells were refused by
the store-seam bands, and the certified route-wave record reads the
block ground as an L=9 state -- 16 pinned components ~90 plus one
71.999 -- the declared 5L is NOT this block's ground; the meters
here re-measure exactly that).

THE SETUP: the trusted ROHF rung loaded; the core (marked: the
certified rotation-closed span); the window = the visible AO columns
of 10f + 13d PLUS the distilled 7s radial -- that one radial is
genuinely distilled production machinery (the 7s shell spans three
s bands and the radial is their occupation-weighted combination), so
it enters through one marked call (``cells.window_columns`` on the
declared key), disclosed here and in the table.  SA over 17
components; the two-pass readout; the engine's own NEVPT per root;
atomterm's meters.

Run from the repo root:
    OMP_NUM_THREADS=16 python scripts/u14e/u14e_vv_6_13.py \
        [--ground-only]  |  tee logs/u14e_vv_6_13.log
"""

import sys
import time
from datetime import date

import numpy as np
from pyscf import mcscf, mrpt

from _vv_common import Verdict, drop_row, load_reference, meters

SECTOR = (5, 1)   # (na, nb): the 7s2 shell's pair joins the 4 opens
NCORE = 5
KEY_613 = "u14e_10f13d7s_5-1_sa-casscf_ao_-_pur_avg"


def window(mol, mf, c_core):
    """10f + 13d visibly (the band-pure AO columns, core-projected,
    Loewdin) + the distilled 7s radial (one marked production call --
    the only production construction in this leg)."""
    from atomterm.cells import CellKey, key_core, window_columns

    labels = [lab.split()[-1] for lab in mol.ao_labels()]
    S = mf.get_ovlp()
    nao = mol.nao
    # the marked call: the distilled span (the 10f/13d bands are the
    # identity; the 7s radial is the production distillation across
    # its three s bands)
    from atomterm.cells import CellKey, key_core, window_columns
    from atomterm.cases import campaign_spec
    spec = campaign_spec("U14e")
    key = CellKey.from_address(KEY_613)
    _nc, _cc, c_pur = key_core(key, mol, mf, spec.shell_bands)
    act, ncas, _rep = window_columns(key, mol, mf, spec.shell_bands,
                                     spec.shell_radials)
    lean = np.max(np.abs(act.T @ S @ c_pur))
    # visible re-orthogonalization against the pur core span (the
    # distilled columns arrive orthogonalized against the canonical
    # prefix -- the same treatment every production cell applies)
    w = act - c_core @ (c_core.T @ S @ act)
    g = w.T @ S @ w
    gv, U = np.linalg.eigh(g)
    if gv.min() < 1e-8:
        sys.exit(f"the distilled window carries a core direction "
                 f"(min gram {gv.min():.2e})")
    c_win = w @ U @ np.diag(gv ** -0.5) @ U.T
    C = np.hstack([c_core, c_win])
    P_all = C @ C.T @ S
    rest = mf.mo_coeff - P_all @ mf.mo_coeff
    g2 = rest.T @ S @ rest
    gv2, U2 = np.linalg.eigh(g2)
    keep = gv2 > 1e-7 * gv2.max()
    c_virt = rest @ U2[:, keep] @ np.diag(gv2[keep] ** -0.5)
    if c_virt.shape[1] != nao - C.shape[1]:
        sys.exit(f"the complement came out {c_virt.shape[1]} wide, "
                 f"expected {nao - C.shape[1]}")
    mo = np.hstack([C, c_virt])
    print(f"  window: the distilled span, {ncas} columns (the "
          f"10f/13d bands are the identity; the 7s radial is the "
          f"production distillation, core-lean {lean:.2f} "
          f"re-orthogonalized away); the window gram's min "
          f"{gv.min():.6f}", flush=True)
    return mo, ncas, lean, float(gv.min())


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    ground_only = "--ground-only" in argv
    v = Verdict()
    mol, mf = load_reference()
    from atomterm.core import purified_core_orbitals
    c_core = purified_core_orbitals(mol, mf, NCORE)[0]
    mo_full, ncas, lean, gv_min = window(mol, mf, c_core)

    t0 = time.time()
    mc = mcscf.CASSCF(mf, ncas, SECTOR)
    mc.state_average_([1.0 / 17] * 17)
    mc.kernel(mo_full)
    print(f"  [sa] ensemble {mc.e_tot:.8f}, converged "
          f"{mc.converged} ({time.time() - t0:.1f} s)", flush=True)
    cas = mcscf.CASCI(mf, ncas, SECTOR)
    cas.mo_coeff = mc.mo_coeff
    cas.canonicalization = False
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
    expect = meters(mol, cas.mo_coeff[:, NCORE:NCORE + ncas], *SECTOR)
    e_tot_0 = float(np.atleast_1d(np.asarray(cas.e_tot)).ravel()[0])
    e_cas = np.atleast_1d(np.asarray(cas.e_cas)).ravel()
    core_e = e_tot_0 - float(e_cas[0])
    rows = []
    for k in (range(1) if ground_only else range(17)):
        pt = mrpt.NEVPT(cas, root=k)
        pt.ci = list(pt.ci)
        e2 = float(pt.kernel())
        rows.append({"root": k, "e_ref": float(e_cas[k] + core_e),
                     "e2": e2, "l2": expect("l2", cas.ci[k]),
                     "s2": expect("s2", cas.ci[k])})
        print(f"  root {k}: e_ref {rows[-1]['e_ref']:.8f}, E(2) "
              f"{e2 * 1e3:+.4f} mEh, <L2> {rows[-1]['l2']:.6f}, "
              f"<S2> {rows[-1]['s2']:.6f}", flush=True)
    r0 = rows[0]
    v.note("the (6,13) meters (the disclosure leg)",
           f"<L2> in [{min(r['l2'] for r in rows):.6f}, "
           f"{max(r['l2'] for r in rows):.6f}], <S2> up to "
           f"{max(r['s2'] for r in rows):.3f} over {len(rows)} "
           f"roots; ensemble {mc.e_tot:.8f} (draw-dependent on this "
           f"flat near-degenerate manifold).  THE FINDING this leg "
           f"exposes: the production-distilled 7s radial is "
           f"{lean:.2f} core-lean and its pur remainder is "
           f"{gv_min:.4f} of the window gram -- the (6,13) "
           f"window's 7s direction is compromised on the pur "
           f"substrate, the class behind the stored (6,13) SA "
           f"refusals (the declared-out record).  The certified "
           f"block-ground classification (L=9) lives in the "
           f"route-wave record's pinned components.  E(2) root 0 "
           f"{r0['e2'] * 1e3:+.4f} mEh, total "
           f"{r0['e_ref'] + r0['e2']:.8f}")
    drop_row({
        "leg": "6_13", "route": "free", "cas": "CAS(6,13)",
        "e_ref": r0["e_ref"], "e2_meh": r0["e2"] * 1e3,
        "e_total": r0["e_ref"] + r0["e2"],
        "anchor_e2_meh": None, "delta_meh": 0.0,
        "l2": f"{min(r['l2'] for r in rows):.6f}.."
              f"{max(r['l2'] for r in rows):.6f}",
        "s2": "L=9 disclosure", "date": date.today().isoformat(),
    })
    v.report()   # notes only -- the disclosure leg never gates
    return 0


if __name__ == "__main__":
    sys.exit(main())
