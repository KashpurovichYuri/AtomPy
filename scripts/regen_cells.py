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
    """The declared key-set: the spine (README · Methods) + the FULL
    manifest grid -- method x mask x substrate x start, the NOs start
    under the radial mask VOID (the seam rejects it: NO columns carry
    no member labels to rotate) and the MO start on pur VOID (the
    start carries its own canonical core; the substrate axis cannot
    enter the run) -- + the virtmix diagnostic and the CASCI pur row.
    The NOs start is the full-span group prefix
    (<window>no2 for {3d,4s}: the d and s occupation groups); k is
    window identity. Deduped by address; P4 folds the historic rung
    names away entirely."""
    sec = tuple(sector)
    keys = spine_keys(atom, window, sec)
    for m in ("casscf", "sa"):
        for mask in (None, "radial", "frozen"):
            for sub in ("pol", "pur"):
                for st in ("ao", "mo", "no"):
                    if mask == "radial" and st != "ao":
                        continue  # radial rotates member-labeled AO slots
                    if st == "mo" and sub == "pur":
                        continue  # VOID at the key: the mo start carries
                        # its own canonical core (substrate moot)
                    keys.append(CellKey(atom, window if st != "no"
                                        else f"{window}no2", sec, method=m,
                                        start=st, mask=mask, substrate=sub))
    keys += [
        CellKey(atom, window, sec, method="casscf", mask="virtmix"),
        CellKey(atom, window, sec, substrate="pur"),  # casci pur row
    ]
    seen: dict[str, CellKey] = {}
    for k in keys:
        seen.setdefault(k.address, k)  # the spine overlaps the grid
    return list(seen.values())


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
    queue = cell_set(name, window, sector)
    if name == "Ti":  # the d² lattice rides the same reference + store
        queue += lattice_keys(name)
    for key in queue:
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
