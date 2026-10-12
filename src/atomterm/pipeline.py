"""End-to-end pipeline: ROHF -> CASCI -> CASSCF -> SA-CASSCF -> L-analysis.

Every stage is checkpointed (see atomterm.drivers); the L-multiplicity of the
obtained states is analyzed with exact Lz / L^2 / S^2 matrices built in the
determinant basis, and symmetry preservation is quantified by ||[H, Lz]||_F
and ||[H, L^2]||_F in each active space.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from atomterm.system import AtomConfig
from atomterm.drivers import (
    build_full_mo,
    run_casci,
    run_casscf,
    run_rohf,
    run_sa_casscf,
    select_shell_ao,
)
from atomterm.operators import build_operators
from atomterm.analysis import analyze, format_report, save_report

__all__ = ["run_carbon", "main"]


def _ncore(cfg: AtomConfig, nelecas: tuple[int, int]) -> int:
    nelectron = cfg.charge
    from pyscf.data import elements

    nelectron += elements.ELEMENTS_PROTON[cfg.symbol.upper().capitalize()]
    return (nelectron - (nelecas[0] + nelecas[1])) // 2


def run_carbon(cfg: AtomConfig | None = None, workdir: str = "results") -> dict:
    cfg = cfg or AtomConfig()
    out = Path(workdir)
    out.mkdir(parents=True, exist_ok=True)
    label = cfg.label

    # ---- stage 1: ROHF (checkpointed) -------------------------------------
    rohf_chk = str(out / f"{label}_rohf.chk")
    mol, mf = run_rohf(cfg, rohf_chk)
    print(f"[ROHF      ] E = {mf.e_tot:.8f} Eh   (chk: {rohf_chk})")

    act_idx = select_shell_ao(mol, "2p")
    ncas = len(act_idx)
    assert ncas == 3, f"expected a 3-orbital 2p shell, got {ncas}"
    iao = np.eye(mol.nao_nr())
    c_shell = iao[:, act_idx]

    results: dict[str, dict] = {"rohf_energy": float(mf.e_tot)}
    reports: list[str] = []

    def stage_analysis(tag: str, cas, c_act, na: int, nb: int):
        ncore = _ncore(cfg, (na, nb))
        ops = build_operators(mol, mf, c_act, ncore, na, nb)
        sector = f"Ms={'+' if na - nb > 0 else ''}{na - nb} (na={na}, nb={nb})"
        res = analyze(ops["H"], ops["Lz"], ops["L2"], ops["S2"], sector)
        e_mine = np.linalg.eigvalsh(ops["H"]) + cas.core_energy
        e_pyscf = np.sort(np.atleast_1d(np.asarray(cas.energies)))[: len(e_mine)]
        dev = float(np.max(np.abs(e_mine[: len(e_pyscf)] - e_pyscf)))
        title = f"[{tag}]  total energies vs PySCF: max|dE| = {dev:.2e} Eh"
        text = format_report(res, title)
        print(text)
        reports.append(text)
        results[tag] = {
            "analysis": res,
            "energies_pyscf": np.asarray(cas.energies).tolist(),
        }

    # ---- stage 2: CASCI on the ROHF solution, triplet sector (Ms=1) -------
    casci1_chk = str(out / f"{label}_casci_ms1.chk")
    ncore1 = _ncore(cfg, (2, 0))
    mo1 = build_full_mo(mf, list(range(ncore1)), act_idx)
    cas1 = run_casci(mf, ncas, (2, 0), mo1, casci1_chk, nroots=3, verbose=0)
    print(f"[CASCI Ms=1] roots: {np.round(cas1.energies, 6)}")
    stage_analysis("casci_ms1", cas1, c_shell, 2, 0)

    # ---- stage 3: CASCI singlet sector (Ms=0): full 1S/1D/3P multiplet ----
    casci0_chk = str(out / f"{label}_casci_ms0.chk")
    ncore0 = _ncore(cfg, (1, 1))
    mo0 = build_full_mo(mf, list(range(ncore0)), act_idx)
    cas0 = run_casci(mf, ncas, (1, 1), mo0, casci0_chk, nroots=9, verbose=0)
    print(f"[CASCI Ms=0] roots: {np.round(cas0.energies, 6)}")
    stage_analysis("casci_ms0", cas0, c_shell, 1, 1)

    # ---- stage 4: CASSCF (single root) on the ROHF solution ---------------
    casscf_chk = str(out / f"{label}_casscf.chk")
    cass = run_casscf(mf, ncas, (2, 0), mo1, casscf_chk, verbose=0)
    print(f"[CASSCF    ] E = {cass.energies[0]:.8f} Eh")
    c_cass = cass.mo_coeff[:, ncore1 : ncore1 + ncas]
    stage_analysis("casscf", cass, c_cass, 2, 0)

    # ---- stage 5: SA-CASSCF over the three 3P M_L components (Ms=1) -------
    sa_chk = str(out / f"{label}_sacasscf.chk")
    sa = run_sa_casscf(mf, ncas, (2, 0), mo1, [1 / 3] * 3, sa_chk, verbose=0)
    print(f"[SA-CASSCF ] roots: {np.round(sa.energies, 6)}")
    c_sa = sa.mo_coeff[:, ncore1 : ncore1 + ncas]
    stage_analysis("sa_casscf", sa, c_sa, 2, 0)

    # ---- persist -----------------------------------------------------------
    (out / f"{label}_report.md").write_text(
        f"# Atomic term analysis: {label}\n\n"
        f"ROHF/cc-pVDZ energy: {mf.e_tot:.8f} Eh\n\n" + "\n\n".join(reports),
        encoding="utf-8",
    )
    import json

    serial = {
        tag: {
            "energies_pyscf": payload["energies_pyscf"],
            "states": [vars(s) for s in payload["analysis"].states],
            "lz_spectrum": payload["analysis"].lz_spectrum,
            "comm_norm_h_lz": payload["analysis"].comm_norm_h_lz,
            "comm_norm_h_l2": payload["analysis"].comm_norm_h_l2,
        }
        for tag, payload in results.items()
        if tag != "rohf_energy"
    }
    (out / f"{label}_report.json").write_text(
        json.dumps({"rohf_energy": results["rohf_energy"], **serial}, indent=2),
        encoding="utf-8",
    )
    print(f"\nreport written to {out}/{label}_report.md and .json")
    return results


def main():
    run_carbon()


if __name__ == "__main__":
    main()
