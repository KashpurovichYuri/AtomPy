"""Regenerate the campaign cells under the CellKey seam: the declared
key-set per atom (``cells.grid_keys`` -- derived from ``cases.CAMPAIGN``,
the one declaration) through ``run_cell``, payloads keyed by cell
address, artifacts ``elements/<Atom>/cells.{json,md}``. Chkfiles land
under ``elements/<Atom>/rungs/<address>.chk`` (gitignored).

Usage:
    python scripts/regen_cells.py [atoms] [selection]

Atoms (default: every declared atom, CAMPAIGN order):
    Ti                      only titanium

Group (--group): one named declared group instead of the full grid
    grid                    the full factor grid, every declared window
                            (the default; plus the declared lattice)
    spine                   the fold spine: raw/avg/frac on the primary
                            window
    diagnostics             the virtual-mix rows only
    planned                 the planned windows' grids (declared,
                            not yet run)
    lattice                 the free-ion d2 lattice cells
    lane                    the pinned-lane legs (CI-level term
                            energies, run_lane_cell on a D2h substrate)

Run gate (--full materializes): cells voided by a declared scenario
    (cases.Scenario -- e.g. the anchor-first protocol on planned
    windows' optimizer cells) are recorded as void-with-reason and
    skipped; --full runs them, re-testing the scenario's dated claim.

Axis filters (AND-composed; each matches the CellKey field exactly):
    --method casci          casci | casscf | sa-casscf
    --window 3d4s           window label, e.g. 2p, 3d4s, 3dno1
    --mask radial           radial | frozen | virtmix | - (maskless)
    --start ao              ao | mo | no
    --substrate pur         pol (canonical core) | pur (closed span)
    --fold avg              raw | avg | frac
    --sector 2-1            the (na, nb) active-electron sector

Store semantics:
    unfiltered run          REPLACES the store (also the pruning path)
    filtered run            MERGES into the existing store (untouched
                            cells keep their stored payloads; a merge
                            never prunes, and a filtered run with no
                            store to merge into is refused -- run
                            unfiltered first)
    the store writes PER-CELL -- a kill mid-group keeps every
    finished cell; cells.md (the render) refreshes at group end.
Threads: OPENBLAS_NUM_THREADS=64 OMP_NUM_THREADS=64 PYSCF_NUM_THREADS=64.
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from pyscf import gto, scf  # noqa: E402

from _tables_common import atomic_write_text  # noqa: E402
from atomterm.cases import CAMPAIGN, campaign_spec, scenarios_for  # noqa: E402
from atomterm.cells import (CellKey, cell_payload, diagnostics,  # noqa: E402
                            grid_keys, lane_keys, lattice_keys, planned,
                            kernel_cell_payload, route_keys, run_cell,
                           run_lane_cell, spine_keys,
                            void_reason)
from atomterm.reference import run_frac_occ, run_rohf  # noqa: E402
from atomterm.system import sector_na_nb  # noqa: E402
from atomterm.windows import select_valence_ao  # noqa: E402


def lane_substrate(spec):
    """The pinned lane's D2h substrate: plain ROHF at tight tolerance
    (the certified lane construction's substrate; the reference
    spec's occupation forcing is deliberately NOT applied -- the lane
    record was measured on the plain substrate)."""
    mol = gto.M(atom=f"{spec.symbol} 0 0 0", basis=spec.basis,
                spin=spec.spin, ecp=spec.ecp, symmetry="D2h", verbose=0)
    mf = scf.ROHF(mol)
    mf.conv_tol = 1e-12
    mf.kernel()
    if not mf.converged:
        sys.exit(f"the lane substrate SCF did not converge "
                 f"({spec.symbol}, D2h)")
    print(f"lane substrate D2h ROHF E = {mf.e_tot:.8f} Eh", flush=True)
    return mol, mf


def run_queue(name, spec, sector, group) -> list[CellKey]:
    """The run's queue: the named declared group, or the full grid
    plus the declared lattice."""
    if group == "spine":
        return spine_keys(name, "".join(spec.shell_patterns), sector)
    if group == "diagnostics":
        return diagnostics(spec)
    if group == "planned":
        return planned(spec)
    if group == "lattice":
        return lattice_keys(spec)
    if group == "lane":
        return lane_keys(spec)
    if group == "route":
        return route_keys(spec)
    return grid_keys(spec) + lattice_keys(spec)


def keep(key: CellKey, args) -> bool:
    """Key-identity filter: every given axis must match (``-`` = the
    None mask, the address grammar's convention)."""
    if args.method and key.method != args.method:
        return False
    if args.window and key.window != args.window:
        return False
    if args.mask and (key.mask or "-") != args.mask:
        return False
    if args.start and key.start != args.start:
        return False
    if args.substrate and key.substrate != args.substrate:
        return False
    if args.fold and key.fold != args.fold:
        return False
    if args.sector and key.sector != args.sector:
        return False
    return True


def parse_args(argv) -> argparse.Namespace:
    """The runner's CLI: atoms, one optional group, axis filters (every
    help string carries its value inventory)."""
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("atoms", nargs="*",
                   help="declared atoms to run (default: every CAMPAIGN "
                        "atom, CAMPAIGN order)")
    p.add_argument("--group", choices=["grid", "spine", "diagnostics",
                                       "planned", "lattice", "lane",
                                       "route"],
                   help="one named declared group instead of the full "
                        "grid: grid | spine | diagnostics | planned | "
                        "lattice | lane | route")
    p.add_argument("--full", action="store_true",
                   help="materialize scenario-voided cells (re-testing "
                        "their dated claims)")
    p.add_argument("--fresh", action="store_true",
                   help="re-run stored cells (ignore the store cache -- "
                        "payload migrations and re-measures)")
    p.add_argument("--method",
                   help="keep this method: casci | casscf | sa-casscf")
    p.add_argument("--window",
                   help="keep this window label, e.g. 2p, 3d4s, 3dno1")
    p.add_argument("--mask",
                   help="keep this mask: radial | frozen | virtmix | - "
                        "(the maskless cells)")
    p.add_argument("--start", help="keep this start: ao | mo | no")
    p.add_argument("--substrate",
                   help="keep this substrate: pol | pur")
    p.add_argument("--fold", help="keep this fold: raw | avg | frac")
    p.add_argument("--sector", metavar="NA-NB",
                   help="keep this (na, nb) sector, e.g. 2-1")
    args = p.parse_args(argv)
    args.sector = (tuple(int(x) for x in args.sector.split("-"))
                   if args.sector else None)
    return args


def render(atom, rohf_e, cells: dict, ruler: str) -> str:
    """The cells.md report: one row per stored FACTS cell
    (key-addressed, the ruler's ΔE within its sector only), plus
    per-root configuration lines; the other payload shapes render in
    their own sections (lane cells: term energies; void cells: their
    scenario's reason).  Markdown pipe tables (the meter names are
    the ladder's; a value column may hold ``nan`` -- no common zero
    for ΔE, or the meter not carried by that payload)."""
    facts_cells = {a: r for a, r in cells.items() if "e_ground" in r}
    other = [(a, r) for a, r in cells.items() if "e_ground" not in r]
    ref = facts_cells.get(ruler)
    head = ("| cell | E_ground | dE_vs_avg | ⟨L²⟩ | ‖[H,Lz]‖ | ‖[H,L²]‖ "
            "| ‖[P,Lz]‖ | ‖[P,L⊥]‖ | win_l | aniso | max dE | ortho | conv |")
    lines = [f"# Orbital-treatment cells: {atom}", "",
             f"ROHF reference energy: {rohf_e:.8f} Eh", "",
             head,
             "|---" * (head.count("|") - 1) + "|"]
    ruler_sec = CellKey.from_address(ruler).sector
    for addr, r in facts_cells.items():
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
            f"| {addr} | {r['e_ground']:.8f} | {d_e:.3e} "
            f"| {r['ground_l2']:.5f} "
            f"| {r['comm_norm_h_lz']:.2e} | {r['comm_norm_h_l2']:.2e} "
            f"| {r['win_lz']:.1e} | {r['win_lperp']:.1e} "
            f"| {purity:.3f} "
            f"| {aniso if aniso is not None else float('nan'):.1e} "
            f"| {r['max_de_vs_pyscf']:.1e} | {r['ortho_err']:.1e} "
            f"| {flag.strip()} |")
    lines += ["", "configurations per root (top-2 determinants):", ""]
    for addr, r in facts_cells.items():
        for i, cfg_str in enumerate(r["configurations"][:9]):
            lines.append(f"- `{addr}` root {i}: {cfg_str}")
    if other:
        lines += ["", "lane cells (tagged windows; the pinned rows "
                      "carry term energies) and void cells "
                      "(scenario-gated):", ""]
        for addr, r in other:
            if "void_reason" in r:
                lines.append(f"- `{addr}` VOID -- {r['void_reason']}")
            elif "total" in r:
                lines.append(
                    f"- `{addr}` E_term {r['total']:.8f}"
                    f"  ptp {r['ptp']:.1e}"
                    f"  L2 {['%.3f' % v for v in r['l2']]}")
            else:                      # a free-route twin (cell payload)
                e = r["energies_pyscf"][0]
                lines.append(f"- `{addr}` E {e:.8f}")
    return "\n".join(lines) + "\n"


def main() -> None:
    args = parse_args(sys.argv[1:])
    unknown = [a for a in args.atoms if a not in CAMPAIGN]
    if unknown:
        sys.exit(f"not declared: {', '.join(unknown)}; "
                 f"declared: {', '.join(CAMPAIGN)}")
    filtered = bool(args.group or args.method or args.window or args.mask
                    or args.start or args.substrate or args.fold
                    or args.sector)
    for name in CAMPAIGN:
        if args.atoms and name not in args.atoms:
            continue
        spec = campaign_spec(name)
        print(f"\n===== cells {name} =====", flush=True)
        out = Path("elements") / name
        (out / "rungs").mkdir(parents=True, exist_ok=True)
        cfg = spec.atom_config()
        sector = sector_na_nb(spec.nelec_act, cfg.spin)
        queue = [k for k in run_queue(name, spec, sector, args.group)
                 if keep(k, args)]
        if not queue:
            print("  no declared cell matches the filter -- skipped",
                  flush=True)
            continue
        cells: dict[str, dict] = {}
        if filtered:
            # merge mode: untouched cells keep their stored payloads.
            # Guards BEFORE any SCF: no store to merge into (a filtered
            # first run would create a subset-only store) and an
            # unreadable store both refuse with the rebuild hint.
            store = out / "cells.json"
            if not store.exists():
                sys.exit(f"no {store} to merge into -- a filtered first "
                         "run would create a subset-only store; run "
                         "unfiltered first")
            try:
                cells = json.loads(store.read_text())["cells"]
            except (json.JSONDecodeError, KeyError) as e:
                sys.exit(f"{store} is unreadable ({e!r}) -- rerun "
                         "unfiltered to rebuild the store")
        mol, mf = run_rohf(cfg, str(out / "rungs" / "rohf.chk"),
                           ref=spec.ref)
        print(f"ROHF E = {mf.e_tot:.8f} Eh", flush=True)
        base = CellKey(name, "".join(spec.shell_patterns), sector)
        act_idx = select_valence_ao(mol, base.shell_patterns)
        _, mf_f = run_frac_occ(cfg, str(out / "rungs" / "fracocc.chk"),
                               spec.nelec_act, len(act_idx),
                               open_idx=act_idx)

        scens = scenarios_for(spec)
        lane_sub = None
        for key in queue:
            reason = void_reason(key, scens)
            if reason and not args.full:
                print(f"  [void] {key.address} -- {reason}", flush=True)
                cells[key.address] = {"cell": key.address,
                                      "void_reason": reason}
                continue
            print(f"  [cell] {key.address}", flush=True)
            if key.fold_tag is not None and key.substrate == "pur":
                if lane_sub is None:
                    lane_sub = lane_substrate(spec)
                cells[key.address] = run_lane_cell(
                    key, lane_sub[0], lane_sub[1], spec.term,
                    fresh=args.fresh, workdir="elements")
            else:
                facts = run_cell(key, mf, frac_mf=mf_f, fresh=args.fresh,
                                 workdir="elements",
                                 budgets={"sa_weights": spec.sa_weights})
                cells[key.address] = (
                    cell_payload(key, facts) if facts is not None
                    else kernel_cell_payload(key, mol, mf))
            # the store writes PER-CELL: a kill mid-group (the shell
            # host caps background jobs) keeps every finished cell
            atomic_write_text(out / "cells.json",
                          json.dumps({"atom": name,
                                      "rohf_energy": float(mf.e_tot),
                                      "cells": cells}, indent=2))
        atomic_write_text(out / "cells.md",
                      render(name, mf.e_tot, cells,
                             ruler=CellKey(name, base.window, sector,
                                           fold="avg").address))
        print(f"  wrote {out}/cells.json + cells.md", flush=True)

    print("\ndone")


if __name__ == "__main__":
    main()
