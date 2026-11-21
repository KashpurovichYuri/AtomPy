"""The reference-comparison V&V leg -- ROHF vs UHF, plus the
spin-broken-start SA question.  THE NEW-PHYSICS LEG: everything is
printed/reported; nothing gates.

  (e) the UHF sibling's spin honesty: E beside the ROHF rung, <S^2>
      by the OVERLAP FORMULA (written out below -- six lines, the
      analysis of interest), the corresponding-orbital pairing
      spectrum;
  the UHF-seeded SA: the visible free-route construction riding the
  UHF object (PySCF's CASSCF accepts it directly) -- does the SA
  repair recover <L2>/<S2> from the spin-broken reference?  Measured
  and reported.

The UHF rung: loaded from ``rungs/uhf.chk`` when present (trusted
cache), else warmed through the production runner (``run_uhf`` with
the basin gate -- one marked call, disclosed).

Run from the repo root:
    OMP_NUM_THREADS=16 python scripts/u14e/u14e_vv_uhf.py \
        |  tee logs/u14e_vv_uhf.log
"""

import sys
import time
from datetime import date
from pathlib import Path

import numpy as np
from pyscf import mcscf, scf

from _vv_common import (ATOM_DIR, BAND_L2, EH_EV, L_L2, S_S2,
                        Verdict, drop_row, load_reference, load_rung,
                        meters)

SECTOR = (4, 0)
NCORE = 5
SHELLS = (("10f", 7), ("13d", 5))


def uhf_rung():
    """The UHF sibling: the trusted chk, else the production runner
    (marked call -- the forced-occupation hooks + the basin gate are
    the reference layer's discipline, not under test here)."""
    chk = str(ATOM_DIR / "rungs" / "uhf.chk")
    if Path(chk).exists():
        mol, mf = load_rung(chk)
        print(f"  UHF (trusted rung, loaded): E = {mf.e_tot:.8f}",
              flush=True)
        return mol, mf
    from atomterm.reference import run_uhf
    from atomterm.cases import campaign_spec
    spec = campaign_spec("U14e")
    mol, mf = run_uhf(spec.atom_config(), chk, ref=spec.ref)
    print(f"  UHF (production runner, marked): E = {mf.e_tot:.8f}",
          flush=True)
    return mol, mf


def contamination(mf):
    """<S^2> by the overlap formula -- written out: Sz(Sz+1) + N_beta
    - sum_ij |<phi_a_i|phi_b_j>|^2 over the occupied blocks; a
    spin-paired determinant makes the sum N_beta and <S^2> = S(S+1)
    exactly; any excess is spin contamination.  Returns
    (s2, sigma) -- sigma the occupied-overlap singular values
    (1 = a paired couple, 0 = fully broken)."""
    ca, cb = mf.mo_coeff
    occ = np.asarray(mf.mo_occ)
    a_idx = np.flatnonzero(occ[0] > 0)
    b_idx = np.flatnonzero(occ[1] > 0)
    na, nb = len(a_idx), len(b_idx)
    x = ca[:, a_idx].conj().T @ mf.get_ovlp() @ cb[:, b_idx]
    sz = (na - nb) / 2.0
    s2 = sz * (sz + 1.0) + nb - float(np.sum(np.abs(x) ** 2))
    return s2, np.linalg.svd(x, compute_uv=False)


def visible_window(mol, mf_rohf, c_core):
    """The same visible window as the 4_12 leg (band-pure shells)."""
    labels = [lab.split()[-1] for lab in mol.ao_labels()]
    sel = [i for i, lab in enumerate(labels)
           if any(lab.startswith(w) for w, _ in SHELLS)]
    if len(sel) != 12:
        sys.exit(f"the label selection found {len(sel)} AOs, "
                 "expected 12")
    nao = mol.nao
    S = mf_rohf.get_ovlp()
    act = np.zeros((nao, len(sel)))
    for k, i in enumerate(sel):
        act[i, k] = 1.0
    P_core = c_core @ c_core.T @ S   # the metric projector
    w = act - P_core @ act
    g = w.T @ S @ w
    gv, U = np.linalg.eigh(g)
    c_win = w @ U @ np.diag(gv ** -0.5) @ U.T
    C = np.hstack([c_core, c_win])
    P_all = C @ C.T @ S   # (C^T S C = I)
    rest = mf_rohf.mo_coeff - P_all @ mf_rohf.mo_coeff
    g2 = rest.T @ S @ rest
    gv2, U2 = np.linalg.eigh(g2)
    keep = gv2 > 1e-7 * gv2.max()
    c_virt = rest @ U2[:, keep] @ np.diag(gv2[keep] ** -0.5)
    if c_virt.shape[1] != nao - C.shape[1]:
        sys.exit(f"the complement came out {c_virt.shape[1]} wide, "
                 f"expected {nao - C.shape[1]}")
    return np.hstack([C, c_virt])


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    v = Verdict()
    mol, mf_rohf = load_reference()
    mol_u, mf_uhf = uhf_rung()
    e_gap = (mf_uhf.e_tot - mf_rohf.e_tot) * 1e3
    print(f"  E_uhf - E_rohf = {e_gap:+.3f} mEh (the broken "
          f"solution's descent)", flush=True)
    s2, sigma = contamination(mf_uhf)
    n_broken = int(np.sum(sigma < 1.0 - 1e-6))
    v.note("(e) the UHF sibling's spin honesty",
           f"<S2> = {s2:.6f} (Sz exact by construction; the exact "
           f"5-fold would read {S_S2:.3f}); pairing spectrum: "
           f"{len(sigma)} occupied pairs, {n_broken} below "
           f"1 - 1e-6, smallest {np.sort(sigma)[:3]}")
    v.note("(e) the ROHF-vs-UHF reference comparison",
           f"ROHF {mf_rohf.e_tot:.8f} vs UHF {mf_uhf.e_tot:.8f} "
           f"(UHF {e_gap:+.3f} mEh lower); the cores assume one "
           f"closed-shell density -- the UHF object is a "
           f"reference-class product, not a cells axis")

    # the UHF-seeded SA: the visible free-route construction riding
    # the UHF object -- does the SA repair recover the term identity
    # from the spin-broken reference?
    from atomterm.core import purified_core_orbitals
    c_core = purified_core_orbitals(mol, mf_rohf, NCORE)[0]
    mo_full = visible_window(mol, mf_rohf, c_core)
    t0 = time.time()
    mc = mcscf.CASSCF(mf_uhf, 12, SECTOR)
    mc.state_average_([1.0 / 17] * 17)
    mc.kernel(mo_full)
    print(f"  [sa] UHF-seeded ensemble {mc.e_tot:.8f}, converged "
          f"{mc.converged} ({time.time() - t0:.1f} s)", flush=True)
    cas = mcscf.CASCI(mf_uhf, 12, SECTOR)
    cas.mo_coeff = mc.mo_coeff
    cas.canonicalization = False
    cas.fcisolver.nroots = 17
    cas.kernel()
    expect = meters(mol_u, cas.mo_coeff[:, NCORE:NCORE + 12],
                    *SECTOR)
    l2s = [expect("l2", cas.ci[k]) for k in range(17)]
    s2s = [expect("s2", cas.ci[k]) for k in range(17)]
    dev = max(max(abs(x - L_L2) for x in l2s),
              max(abs(x - S_S2) for x in s2s))
    if dev <= BAND_L2:
        v.note("the UHF-seeded SA recovers <L2> = 72 / <S2> = 6",
               f"{len(l2s)} roots from the spin-broken reference, max "
               f"dev {dev:.2e} (band {BAND_L2:g}) -- the SA repair "
               f"answers the recovery question")
    else:
        v.note("the UHF-seeded SA from the RAW start descends "
               "broken",
               f"<L2> dev {dev:.2e} (the term identity lost; <L2> in "
               f"[{min(l2s):.4f}, {max(l2s):.4f}], ensemble "
               f"{mc.e_tot:.4f} -- an energy of the UHF-CASCI "
               f"partition, not comparable to the ROHF-route optimum "
               f"-50.4238).  The committed orchestrator record "
               f"(logs/u14e_vv.log) measured the CONVERGED-start "
               f"reading: the repair holds there (dev 8.9e-15, 2 "
               f"cycles).  Together: the SA repair holds the "
               f"converged gauge; a raw spin-broken start escapes")
    drop_row({
        "leg": "ref", "route": "uhf", "cas": "CAS(4,12)",
        "e_ref": float(mf_uhf.e_tot), "e2_meh": 0.0,
        "e_total": float(mf_uhf.e_tot),
        "anchor_e2_meh": None, "delta_meh": e_gap,
        "l2": f"{min(l2s):.6f}..{max(l2s):.6f}",
        "s2": f"<S2> {s2:.4f} (UHF det)", "date": date.today().isoformat(),
    })
    v.report()
    return 0


if __name__ == "__main__":
    sys.exit(main())
