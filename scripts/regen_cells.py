"""Regenerate the campaign cells under the CellKey seam: the declared
key-set per atom (``cells.grid_keys`` -- derived from ``cases.CAMPAIGN``,
the one declaration) through ``run_cell``, payloads keyed by cell
address, artifacts ``elements/<Atom>/cells.{json,md}``. Chkfiles land
under ``elements/<Atom>/rungs/<address>.chk`` (gitignored).

Usage: python scripts/regen_cells.py [Ti C Si Sc]
Threads: OPENBLAS_NUM_THREADS=64 OMP_NUM_THREADS=64 PYSCF_NUM_THREADS=64.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, "src")
from atomterm.cases import CAMPAIGN
from atomterm.cells import CellKey, cell_payload, grid_keys, run_cell
from atomterm.reference import run_frac_occ, run_rohf
from atomterm.system import sector_na_nb
from atomterm.windows import select_valence_ao


def lattice_keys(atom) -> list[CellKey]:
    """The d² lattice cells (Ti): the free-ion series CAS(2,·) -- the
    4s² pair closed in the core (ncore=10), sector (2,0) = Ms=1, over
    the growing m-complete shell spans. Constructions per row: the raw
    shell AOs (exactly closed), the PySCF canonical MO slice (the
    m-delocalized control, fact 22d), and the full-span group-prefix
    NOs of the source's own term-averaged density (fact 23: k=K spans
    the whole source window -- NO == AOs by construction, the closure
    certificate is the point; k=1 for {3d} doubles as the machinery
    pin). All sectors are tiny (<= 78 dets): seconds each."""
    sec = (2, 0)
    keys = []
    for w in ("3d", "3d4p", "3d4p4d"):
        keys.append(CellKey(atom, w, sec, method="casci", fold="avg"))
        keys.append(CellKey(atom, w, sec, method="casci", start="mo",
                            fold="avg"))
    for w in ("3dno1", "3d4pno2", "3d4p4dno3"):
        keys.append(CellKey(atom, w, sec, method="casci", start="no",
                            fold="avg"))
    return keys


def render(atom, rohf_e, cells: dict, ruler: str) -> str:
    ref = cells.get(ruler)
    lines = [f"# Orbital-treatment cells: {atom}", "",
             f"ROHF reference energy: {rohf_e:.8f} Eh", "",
             f"{'cell':<44} {'E_ground':>16} {'dE_vs_avg':>10} {'<L^2>':>9} "
             f"{'|[H,Lz]|':>10} {'|[H,L2]|':>10} {'|[P,Lz]|':>9} {'|[P,Lp]|':>9} "
             f"{'win_l':>6} {'aniso':>9} {'max|dE|':>9} {'ortho':>9}  conv"]
    ruler_sec = CellKey.from_address(ruler).sector
    for addr, r in cells.items():
        # ΔE vs the ruler only within its sector -- cross-sector rows
        # (the d² lattice vs the 3d²4s² spine) have no common zero
        d_e = (r["e_ground"] - ref["e_ground"]
               if ref and CellKey.from_address(addr).sector == ruler_sec
               else float("nan"))
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
for name, spec in CAMPAIGN.items():
    if only and name not in only:
        continue
    print(f"\n===== cells {name} =====", flush=True)
    out = Path("elements") / name
    (out / "rungs").mkdir(parents=True, exist_ok=True)
    cfg = spec.atom_config()
    mol, mf = run_rohf(cfg, str(out / "rungs" / "rohf.chk"), ref=spec.ref)
    print(f"ROHF E = {mf.e_tot:.8f} Eh", flush=True)
    sector = sector_na_nb(spec.nelec_act, cfg.spin)
    base = CellKey(name, "".join(spec.shell_patterns), sector)
    act_idx = select_valence_ao(mol, base.shell_patterns)
    _, mf_f = run_frac_occ(cfg, str(out / "rungs" / "fracocc.chk"),
                               spec.nelec_act, len(act_idx),
                               open_idx=act_idx)

    cells: dict[str, dict] = {}
    queue = grid_keys(spec)     # the one declaration (cases.CAMPAIGN)
    if name == "Ti":  # the d² lattice rides the same reference + store
        queue += lattice_keys(name)
    for key in queue:
        print(f"  [cell] {key.address}", flush=True)
        facts = run_cell(key, mf, frac_mf=mf_f, workdir="elements",
                         budgets={"sa_weights": spec.sa_weights})
        cells[key.address] = cell_payload(key, facts)

    (out / "cells.json").write_text(
        json.dumps({"atom": name, "rohf_energy": float(mf.e_tot),
                    "cells": cells}, indent=2), encoding="utf-8")
    (out / "cells.md").write_text(
        render(name, mf.e_tot, cells,
               ruler=CellKey(name, base.window, sector,
                             fold="avg").address),
        encoding="utf-8")
    # the rungs.json/md era ends here: content is key-addressed now
    for ext in ("json", "md"):
        old = out / f"rungs.{ext}"
        if old.exists():
            old.unlink()
    print(f"  wrote {out}/cells.json + cells.md", flush=True)

print("\ndone")
