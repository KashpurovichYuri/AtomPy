"""Ticket 25 diagnostic: our h_eff vs PySCF's h1eff, and the
window-span hypothesis.

Three comparisons on the Ti gated reference, CAS(4,6)/(3,1), shell
window {3d,4s}:
  A. h1(ours, RAW fold, the RUN's columns) vs PySCF's get_h1eff on the
     same columns  -> is there a fold-CONVENTION difference at all?
  B. h1(ours, AVG fold, run columns) vs PySCF's               -> the
     averaging-only difference (expected: the 2.8 mEh class).
  C. the CI in the RAW shell AO span (stage_facts' CASCI-stage
     convention) vs the CI in the RUN's projected columns      -> the
     window-definition offset (the 4s|3s-core overlap hypothesis).
Carbon (closed shell, 2p window orthogonal to core) = the negative
control: everything should agree to ~1e-12.

Usage: python scripts/diag_fold.py
"""

import sys
from pathlib import Path

sys.path.insert(0, "src")
import numpy as np
from pyscf import mcscf

from atomterm.system import AtomConfig
from atomterm import drivers
from atomterm.operators import active_space_integrals
from atomterm.detbasis import determinant_basis, slater_condon_matrix

TI_REF = drivers.ReferenceSpec(
    occ_closed=("1s", "2s", "3s", "4s", "2p", "3p"),
    occ_open=("3d",),
    occ_e_ref=-848.40648875,
)
from run_u import BAS, ECP


def compare(label, mol, mf, act_idx, ncore, na, nb):
    out = Path("results/diag"); out.mkdir(exist_ok=True)
    chk = str(out / f"{label}_diag.chk")
    mo = drivers.build_full_mo(mf, drivers.occupied_core_idx(mf, ncore),
                               act_idx)
    ncas = len(act_idx)
    cas = drivers.run_casci(mf, ncas, (na, nb), mo, chk, nroots=1)
    c_run = cas.mo_coeff[:, ncore:ncore + ncas]
    c_raw = np.eye(mol.nao_nr())[:, act_idx]

    h1_rawfold, g_run = active_space_integrals(
        mol, mf, c_run, ncore, spherical_core=False,
        c_core=cas.mo_coeff[:, :ncore])
    h1_avgfold, _ = active_space_integrals(
        mol, mf, c_run, ncore, spherical_core=True,
        c_core=cas.mo_coeff[:, :ncore])

    mc = mcscf.CASCI(mf, ncas, (na, nb))
    mc.mo_coeff = cas.mo_coeff
    h1p, ecorep = mc.get_h1eff()

    from atomterm.operators import build_operators
    def e0(c_act, spherical):
        ops = build_operators(mol, mf, c_act, ncore, na, nb,
                              spherical_core=spherical,
                              c_core=cas.mo_coeff[:, :ncore])
        return float(np.linalg.eigvalsh(ops["H"].real)[0])

    print(f"\n===== {label} (ncore={ncore}, ncas={ncas}) =====")
    print(f"PySCF e_cas          = {cas.mc.e_cas if hasattr(cas.mc, 'e_cas') else float(cas.mc.e_tot):.9f}")
    print(f"A |h1_ours_raw - h1_pyscf|  max = {np.abs(h1_rawfold - h1p).max():.3e}   (convention test)")
    print(f"B |h1_ours_avg - h1_pyscf|  max = {np.abs(h1_avgfold - h1p).max():.3e}   (averaging test)")
    e_run, e_ours_avg = e0(c_run, False), e0(c_run, True)
    print(f"e_cas(pyscf cols, pyscf h1)  = {e_run:.9f}")
    print(f"e_cas(our build, avg fold)  = {e_ours_avg:.9f}   dE = {(e_ours_avg - e_run) * 1e3:+.3f} mEh  (averaging effect)")

    # C: the same CI in the RAW shell AO span (stage_facts' CASCI
    # convention), averaged ruler
    e_rawspan = e0(c_raw, True)
    print(f"C e_cas(RAW AO span, our avg)= {e_rawspan:.9f}   "
          f"dE vs run = {(e_rawspan - e_run) * 1e3:+.3f} mEh  (window-definition test)")
    ov = float(np.abs(c_run.conj().T @ mf.get_ovlp() @ c_raw
                      - np.eye(ncas)).max())
    s = mf.get_ovlp()
    cross = float(np.abs(c_raw.conj().T @ s @ cas.mo_coeff[:, :ncore]).max())
    print(f"    run-vs-raw span deviation |C_run^T S C_raw - 1| = {ov:.3e};"
          f" max |raw window x core| overlap = {cross:.3f}")


def main():
    Path("results/diag").mkdir(parents=True, exist_ok=True)
    # gated U reference (ticket 21): whole-l spans + Newton
    U_REF = drivers.ReferenceSpec(
        occ_closed=("l0", "l1"), occ_open=("l2", "l3"),
        occ_e_ref=-50.38185233, conv_accelerator="newton")
    U_BANDS = {"f": (0.5, 2.0), "d": (0.1, 0.5)}

    cases = [
        ("C", AtomConfig("C", "ccpvdz", spin=2), None,
         ["2p"], 2, 1, 1, 1),
        ("Si", AtomConfig("Si", "ccpvdz", spin=2), None,
         ["3p"], 2, 1, 1, 6),  # core = [Ne] (1s2s2p), not 1s2s
        ("Sc", AtomConfig("Sc", "ccpvdz", spin=1), None,
         ["3d", "4s"], 3, 2, 1, 9),
        ("Ti", AtomConfig("Ti", "ccpvdz", spin=2), TI_REF,
         ["3d", "4s"], 4, 3, 1, 9),
        ("U", AtomConfig("U", BAS, ecp=ECP, spin=4, tag="U_stuttgart78_dyall"),
         U_REF, None, 4, 4, 0, 5),  # window = 2f+2d bands (built below)
    ]
    for (label, cfg, ref, patterns, nelec, na, nb, ncore) in cases:
        print(f"\n{label} (ROHF{' gated' if ref else ''}) ...", flush=True)
        mol, mf = drivers.run_rohf(cfg, f"results/diag/{label}_rohf.chk",
                                   ref=ref)
        if patterns is not None:
            act = drivers.select_valence_ao(mol, patterns)
        else:
            f_rows = drivers.select_primitives(mol, 3, exponent_band=U_BANDS["f"])
            d_rows = drivers.select_primitives(mol, 2, exponent_band=U_BANDS["d"])
            act = f_rows + d_rows
        compare(label, mol, mf, act, ncore, na, nb)


if __name__ == "__main__":
    main()
