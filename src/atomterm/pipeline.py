"""End-to-end pipeline: ROHF -> CASCI -> CASSCF -> SA-CASSCF -> L-analysis.

Every stage is checkpointed (see atomterm.drivers); the L-multiplicity of the
obtained states is analyzed with exact Lz / L^2 / S^2 matrices built in the
determinant basis, and symmetry preservation is quantified by ||[H, Lz]||_F
and ||[H, L^2]||_F in each active space.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from atomterm.system import AtomConfig
from atomterm.drivers import (
    build_full_mo,
    run_casci,
    run_casscf as casscf_stage,
    run_rohf,
    run_sa_casscf,
    select_valence_ao,
)
from atomterm.stages import stage_facts
from atomterm.analysis import format_report

__all__ = ["run_atom", "run_carbon", "main"]


def _ncore(cfg: AtomConfig, nelecas: tuple[int, int]) -> int:
    nelectron = cfg.charge
    from pyscf.data import elements

    nelectron += elements.ELEMENTS_PROTON[cfg.symbol.upper().capitalize()]
    return (nelectron - (nelecas[0] + nelecas[1])) // 2


def _ndet(ncas: int, na: int, nb: int) -> int:
    """Determinants in a sector: C(norb, na) * C(norb, nb)."""
    return math.comb(ncas, na) * math.comb(ncas, nb)


def _ms_label(na: int, nb: int) -> str:
    """Filename/report label of a sector's M_S = (na - nb) / 2.

    Whole numbers stay bare ("1", "0" -- matching the historic carbon chkfile
    names); halves insert a "p" so paths stay clean ("0p5" for a doublet).
    """
    assert na >= nb, "sectors are built with na >= nb"
    two_ms = na - nb
    if two_ms % 2 == 0:
        return str(two_ms // 2)
    return f"{two_ms / 2:.1f}".replace(".", "p")


def run_atom(
    cfg: AtomConfig,
    shell_patterns: list[str],
    nelec_act: int,
    workdir: str = "results",
    sa_weights: list[float] | None = None,
    run_casscf: bool = True,
    nroots_primary: int = 9,
) -> dict:
    """End-to-end pipeline for any single atom.

    shell_patterns : valence-shell label prefixes selecting the active AO
        space, concatenated in order (e.g. ``["2p"]`` for C 2p^2,
        ``["3d", "4s"]`` for Sc 3d^1 4s^2). Each pattern must match at least
        one AO; each should contribute a complete m-set (rotation closure).
    nelec_act : electrons treated as active. The Ms sectors follow from
        cfg.spin = 2S: the primary sector (na, nb) = ((n+2S)/2, (n-2S)/2),
        plus an Ms=0 sector (n/2, n/2) when n is even and S > 0 (skipped when
        odd -- no Ms=0 exists, e.g. a doublet).
    sa_weights : run SA-CASSCF over these state weights (primary sector).
    run_casscf : run the single-root CASSCF stage (off for transition metals,
        where orbital optimization is multi-solution).
    nroots_primary : PySCF roots for the primary CASCI (clamped to the sector
        size; the Ms=0 CASCI always diagonalizes the full sector so the
        complete multiplet table is produced).

    Every stage is checkpointed under ``workdir`` and analyzed with the exact
    L/S operator matrices; see :func:`run_carbon` for the carbon instance.
    """
    out = Path(workdir)
    out.mkdir(parents=True, exist_ok=True)
    label = cfg.label
    spin = cfg.spin  # 2S (PySCF convention)
    assert spin >= 0 and (nelec_act + spin) % 2 == 0, (
        f"spin 2S={spin} and nelec_act={nelec_act} must have equal parity"
    )

    # ---- stage 1: ROHF (checkpointed) -------------------------------------
    rohf_chk = str(out / f"{label}_rohf.chk")
    mol, mf = run_rohf(cfg, rohf_chk)
    print(f"[ROHF      ] E = {mf.e_tot:.8f} Eh   (chk: {rohf_chk})")

    act_idx = select_valence_ao(mol, shell_patterns)
    ncas = len(act_idx)
    c_shell = np.eye(mol.nao_nr())[:, act_idx]

    # Ms sectors: (na, nb, nroots) for each CASCI diagonalization
    na_p, nb_p = (nelec_act + spin) // 2, (nelec_act - spin) // 2
    sector_defs = [(na_p, nb_p, min(nroots_primary, _ndet(ncas, na_p, nb_p)))]
    if nelec_act % 2 == 0 and spin > 0:
        na0 = nelec_act // 2
        sector_defs.append((na0, na0, _ndet(ncas, na0, na0)))

    results: dict[str, dict] = {"rohf_energy": float(mf.e_tot)}
    reports: list[str] = []

    def stage_analysis(tag: str, cas, c_act, na: int, nb: int):
        ncore = _ncore(cfg, (na, nb))
        # fold with the stage's OWN core (inside stage_facts: it uses
        # cas.mo_coeff[:, :ncore]): orbital optimization may legally
        # reorganize the inactive span, and a different fold than the one
        # PySCF's energies used breaks E_tot = E_core + e_cas. For CASCI
        # stages the passed core IS the mf core (built by build_full_mo).
        sector = f"Ms={'+' if na - nb > 0 else ''}{na - nb} (na={na}, nb={nb})"
        facts = stage_facts(mol, mf, cas, c_act, ncore, na, nb, tag=sector)
        title = (f"[{tag}]  total energies vs PySCF: "
                 f"max|dE| = {facts.max_de_vs_pyscf:.2e} Eh")
        text = format_report(facts.analysis, title)
        off_shell = 1.0 - facts.window_rows.max(axis=1)
        text += (f"\nwindow composition (l-weights): "
                 + ", ".join(f"{n}: {v:.3f}" for n, v in
                             facts.window.items() if v > 1e-4)
                 + f"   (mean off-dominant-l weight {off_shell.mean():.2e};"
                   " l-purity is necessary, not sufficient, for closure)")
        print(text)
        reports.append(text)
        results[tag] = {
            "analysis": facts.analysis,
            "energies_pyscf": np.asarray(facts.energies_pyscf).tolist(),
            "window_weights": facts.window,
            "max_de_vs_pyscf": facts.max_de_vs_pyscf,
            # particle-number sanity of the built determinant space: N is
            # diagonal by construction and must equal na + nb everywhere
            "n_dev": facts.n_dev,
        }

    # ---- stage 2: CASCI in every Ms sector (checkpointed) ------------------
    mo_by_sector: dict[tuple[int, int], np.ndarray] = {}
    for na, nb, nroots in sector_defs:
        ncore = _ncore(cfg, (na, nb))
        mo = build_full_mo(mf, list(range(ncore)), act_idx)
        mo_by_sector[(na, nb)] = mo
        tag = f"casci_ms{_ms_label(na, nb)}"
        cas = run_casci(mf, ncas, (na, nb), mo, str(out / f"{label}_{tag}.chk"),
                        nroots=nroots, verbose=0)
        print(f"[CASCI Ms={_ms_label(na, nb)}] roots: {np.round(cas.energies, 6)}")
        stage_analysis(tag, cas, c_shell, na, nb)

    mo_p = mo_by_sector[(na_p, nb_p)]
    ncore_p = _ncore(cfg, (na_p, nb_p))

    # ---- stage 3: CASSCF (single root) on the primary sector ---------------
    if run_casscf:
        casscf_chk = str(out / f"{label}_casscf.chk")
        cass = casscf_stage(mf, ncas, (na_p, nb_p), mo_p, casscf_chk, verbose=0)
        print(f"[CASSCF    ] E = {cass.energies[0]:.8f} Eh")
        c_cass = cass.mo_coeff[:, ncore_p : ncore_p + ncas]
        stage_analysis("casscf", cass, c_cass, na_p, nb_p)

    # ---- stage 4: SA-CASSCF over the requested state weights ---------------
    if sa_weights:
        sa_chk = str(out / f"{label}_sacasscf.chk")
        sa = run_sa_casscf(mf, ncas, (na_p, nb_p), mo_p, list(sa_weights),
                           sa_chk, verbose=0)
        print(f"[SA-CASSCF ] roots: {np.round(sa.energies, 6)}")
        c_sa = sa.mo_coeff[:, ncore_p : ncore_p + ncas]
        stage_analysis("sa_casscf", sa, c_sa, na_p, nb_p)

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
            "window_weights": payload["window_weights"],
            "max_de_vs_pyscf": payload["max_de_vs_pyscf"],
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


def run_carbon(cfg: AtomConfig | None = None, workdir: str = "results") -> dict:
    """Carbon pilot (2p^2, CAS(2,3)): the generic pipeline with carbon data."""
    return run_atom(
        cfg or AtomConfig(),
        ["2p"],
        2,
        workdir=workdir,
        sa_weights=[1 / 3] * 3,  # SA over the three 3P M_L components (Ms=1)
    )


def main():
    run_carbon()


if __name__ == "__main__":
    main()
