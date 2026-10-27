"""Regenerate the campaign cells under the CellKey seam (the cell wave,
P3): the declared key-set per atom through ``run_cell``, payloads keyed
by cell address, artifacts ``elements/<Atom>/cells.{json,md}`` (the
rungs.json rename -- content became key-addressed). Chkfiles land under
``elements/<Atom>/rungs/<address>.chk`` (gitignored).

Usage: python scripts/regen_cells.py [Ti C Si Sc]
Threads: OPENBLAS_NUM_THREADS=64 OMP_NUM_THREADS=64 PYSCF_NUM_THREADS=64.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, "src")
from atomterm.cells import CellKey, cell_payload, run_cell, spine_keys
from atomterm.reference import ReferenceSpec, run_rohf, run_frac_occ_scf
from atomterm.system import AtomConfig, sector_occupations
from atomterm.windows import select_valence_ao

# the atom's reference data (AGENTS fact 21; Ti identical to tests/conftest.TI_REF;
# Sc wired after a broken-basin fresh draw 2026-10-27: only 5 closed MOs)
REFS = {
    "Ti": ReferenceSpec(occ_closed=("1s", "2s", "3s", "4s", "2p", "3p"),
                        occ_open=("3d",), occ_e_ref=-848.40648875),
    "Sc": ReferenceSpec(occ_closed=("1s", "2s", "2p", "3s", "3p", "4s"),
                        occ_open=("3d",), occ_e_ref=-759.73625510),
}

CASES = {
    "C": (AtomConfig(), "2p", 2, [1 / 3] * 3, 2, 3),
    "Si": (AtomConfig("Si", "ccpvdz", spin=2), "3p", 2, [1 / 3] * 3, 2, 3),
    "Sc": (AtomConfig("Sc", "ccpvdz", spin=1), "3d4s", 3, [1 / 5] * 5, 3, 6),
    "Ti": (AtomConfig("Ti", "ccpvdz", spin=2), "3d4s", 4, [1 / 7] * 7, 4, 6),
}


def cell_set(atom, window, sector) -> list[CellKey]:
    """The declared key-set: the spine (README · Methods) + the
    diagnostics (radial x SA gauge anchor on both substrates, the MO
    start, virtmix, the pur CASCI window-definition row). Mirrors the
    historic rung registry one-to-one; P4 folds the registry away."""
    sec = tuple(sector)
    keys = spine_keys(atom, window, sec)
    keys += [
        CellKey(atom, window, sec, method="sa"),       # sa free, pol
        CellKey(atom, window, sec, method="sa", mask="radial"),
        CellKey(atom, window, sec, method="sa", mask="radial",
                substrate="pur"),
        CellKey(atom, window, sec, method="casscf", start="mo"),
        CellKey(atom, window, sec, method="sa", start="mo"),
        CellKey(atom, window, sec, method="casscf", mask="radial"),
        CellKey(atom, window, sec, method="casscf", mask="virtmix"),
        CellKey(atom, window, sec, method="casscf", mask="frozen"),
        CellKey(atom, window, sec, substrate="pur"),  # casci pur row
        CellKey(atom, window, sec, method="casscf", mask="radial",
                substrate="pur"),
    ]
    seen: dict[str, CellKey] = {}
    for k in keys:
        seen.setdefault(k.address, k)  # spine may overlap diagnostics
    return list(seen.values())


def render(atom, rohf_e, cells: dict, ruler: str) -> str:
    ref = cells.get(ruler)
    lines = [f"# Orbital-treatment cells: {atom}", "",
             f"ROHF reference energy: {rohf_e:.8f} Eh", "",
             f"{'cell':<44} {'E_ground':>16} {'dE_vs_avg':>10} {'<L^2>':>9} "
             f"{'|[H,Lz]|':>10} {'|[H,L2]|':>10} {'|[P,Lz]|':>9} {'|[P,Lp]|':>9} "
             f"{'win_l':>6} {'aniso':>9} {'max|dE|':>9} {'ortho':>9}  conv"]
    for addr, r in cells.items():
        d_e = r["e_ground"] - ref["e_ground"] if ref else 0.0
        purity = max(r["window"].values())
        aniso = r.get("core_anisotropy")
        flag = ""
        if r["diverged"]:
            flag += " DIVERGED"
        if r.get("converged", True) is False:
            flag += " NOT-CONVERGED"
        lines.append(
            f"{addr:<44} {r['e_ground']:>16.8f} {d_e:>10.3e} "
            f"{r['ground_l2']:>9.5f} "
            f"{r['comm_norm_h_lz']:>10.2e} {r['comm_norm_h_l2']:>10.2e} "
            f"{r['win_lz']:>9.1e} {r['win_lperp']:>9.1e} "
            f"{purity:>6.3f} {aniso if aniso is not None else float('nan'):>9.1e} "
            f"{r['max_de_vs_pyscf']:>9.1e} {r['ortho_err']:>9.1e}"
            f"{flag}")
    lines += ["", "configurations per root (top-2 determinants):"]
    for addr, r in cells.items():
        for i, cfg_str in enumerate(r["configurations"][:9]):
            lines.append(f"  {addr:<44} root {i}: {cfg_str}")
    return "\n".join(lines) + "\n"


only = sys.argv[1:] or None
for name, (cfg, window, nelec, saw, fo, no) in CASES.items():
    if only and name not in only:
        continue
    print(f"\n===== cells {name} =====", flush=True)
    out = Path("elements") / name
    (out / "rungs").mkdir(parents=True, exist_ok=True)
    mol, mf = run_rohf(cfg, str(out / "rungs" / "rohf.chk"), ref=REFS.get(name))
    print(f"ROHF E = {mf.e_tot:.8f} Eh", flush=True)
    sector = sector_occupations(nelec, cfg.spin)
    base = CellKey(name, window, sector)
    act_idx = select_valence_ao(mol, base.shell_patterns)
    _, mf_f = run_frac_occ_scf(cfg, str(out / "rungs" / "fracocc.chk"),
                               fo, no, open_idx=act_idx)

    cells: dict[str, dict] = {}
    for key in cell_set(name, window, sector):
        print(f"  [cell] {key.address}", flush=True)
        facts = run_cell(key, mf, frac_mf=mf_f, workdir="elements",
                         budgets={"sa_weights": saw})
        cells[key.address] = cell_payload(key, facts)

    (out / "cells.json").write_text(
        json.dumps({"atom": name, "rohf_energy": float(mf.e_tot),
                    "cells": cells}, indent=2), encoding="utf-8")
    (out / "cells.md").write_text(
        render(name, mf.e_tot, cells,
               ruler=CellKey(name, window, sector, fold="avg").address),
        encoding="utf-8")
    # the rungs.json/md era ends here: content is key-addressed now
    for ext in ("json", "md"):
        old = out / f"rungs.{ext}"
        if old.exists():
            old.unlink()
    print(f"  wrote {out}/cells.json + cells.md", flush=True)

print("\ndone")
