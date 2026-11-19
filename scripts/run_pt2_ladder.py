"""The PT2 ladder collector: NEVPT2 rows over the DECLARED cell set,
written to ``elements/<Atom>/perturbation.json`` -- the renderer's
only source (the render contract: the table and its rows cannot
diverge).

THE DECLARED SET (the intent, one list per atom -- stale cells cannot
leak in; a missing address refuses before any run):
  the span ladder     CASCI {compact, @r1, @s1[, @s2]} x {ao, mo, moc}
  the NO twins        CASCI {@r1no1, @s1no1, 2pno1[, @s2no1]}
  the reference axis  {casscf, sa-casscf} x {compact, @r1, @s1}
                      x {ao, mo}
  the frozen pair     the core-relaxation disclosure rows
The LANE (the core-partition axis) is DECLARED OUT for now: the lane
payload is its own store class (per-irrep-block pinned solutions,
term energies -- no ci_stack/core_energy), so the stored-roots seam
needs a lane-native readout twin before it can price lane rows
(ticket 45's follow-up).
Every stored root is a row (the declared-roots policy -- no silent
draws).  Per-root meters (L2 / Lz / S2) are measured LIVE through the
CI-order bridge; cells.json is never a source.  A stored-roots
refusal does not abort the ladder: the cell records its refusal
(reason string) and the ladder moves on -- the recovery loop is
``scripts/regen_cells.py`` scoped to the refused address, then a
collector re-run fills the row.

Usage:
    python scripts/run_pt2_ladder.py <atom> [<atom> ...]
"""

import json
import sys
import time
from datetime import date, datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402
from pyscf import __version__ as pyscf_version  # noqa: E402

from atomterm.cases import cell_substrate  # noqa: E402
from atomterm.cells import CellKey, cell_chk  # noqa: E402
from atomterm.operators import build_operators, ci_in_detbasis_order  # noqa: E402
from atomterm.perturb import nevpt2_rows  # noqa: E402

EH_EV = 27.211386245988


def _declared_addresses(atom: str) -> list[str]:
    """The declared cell set per atom (the manifest -- edit here)."""
    head, tagged = {
        "C": ("c_2p", ["c_2p@r1", "c_2p@s1"]),
        "Si": ("si_3p", ["si_3p@r1", "si_3p@s1", "si_3p@s2"]),
    }[atom]
    windows = [head] + tagged
    ladder = [f"{x}_2-0_casci_{s}_-_pur_avg"
              for x in windows for s in ("ao", "mo", "moc")]
    no_tails = ([f"{head}@r1no1", f"{head}@s1no1", f"{head}no1"]
                + ([f"{head}@s2no1"] if atom == "Si" else []))
    ladder += [f"{x}_2-0_casci_no_-_pur_avg" for x in no_tails]
    ref_axis = [f"{x}_2-0_{m}_{s}_-_pur_avg"
                for x in windows
                for m in ("casscf", "sa-casscf")
                for s in ("ao", "mo")]
    lane = {"C": "c_2p@1s2s", "Si": "si_3p@1s2s2p3s"}[atom]
    del lane  # the lane axis: DECLARED OUT (its payload class has no
    # ci_stack/core_energy -- the stored-roots seam needs a lane-native
    # readout twin first; ticket 45's follow-up)
    frozen = [f"{x}_2-0_sa-casscf_ao_frozen_pur_avg"
              for x in (head, f"{head}@r1")]
    return ladder + ref_axis + frozen


def _expect(mat, ci_pyscf, norb, na, nb) -> float:
    w = ci_in_detbasis_order(np.asarray(ci_pyscf).ravel(), norb, na, nb)
    return float(np.real(np.vdot(w, mat @ w)))


def collect_atom(atom: str, workdir: str = "elements") -> dict:
    """One atom's declared ladder: rows + meters + provenance."""
    missing = [a for a in _declared_addresses(atom)
               if not cell_chk(CellKey.from_address(a), workdir).exists()]
    if missing:
        sys.exit(f"{atom}: {len(missing)} declared addresses are not in "
                 f"the store: {missing[:3]}{'...' if len(missing) > 3 else ''}")
    entries = []
    for addr in _declared_addresses(atom):
        key = CellKey.from_address(addr)
        na, nb = key.sector
        t0 = time.time()
        mol, mf, data = cell_substrate(key, workdir)
        ncas = int(data.ncas)
        nroots = len(data.ci)
        ncore = (mol.nelectron - (na + nb)) // 2
        mo = np.asarray(data.mo_coeff)
        c_act = mo[:, ncore:ncore + ncas]
        ops = build_operators(mol, mf, c_act, ncore, na, nb,
                              spherical_core=False, c_core=mo[:, :ncore])
        entry = {
            "address": addr,
            "ncas": ncas, "nroots": nroots,
            "sector": [na, nb],
            "core_energy": float(data.core_energy),
            "kernel_energies": np.asarray(
                data.energies, dtype=float).ravel().tolist(),
            "converged": data.converged,
            "cycles": data.cycles,
            "roots": [], "disposition": "rows", "refusal_reason": None,
        }
        try:
            rows = nevpt2_rows(mol, mf, mo, ncas, (na, nb), nroots,
                               ref_energies=data.energies, ref_ci=data.ci)
        except ValueError as exc:
            entry["disposition"] = "refusal"
            entry["refusal_reason"] = str(exc)
            print(f"  [refused] {addr}: {exc}", flush=True)
            entries.append(entry)
            continue
        for r, row in enumerate(rows):
            ci_p = np.asarray(data.ci[r]).ravel()
            entry["roots"].append({
                "root": r,
                "e_ref": row.e_ref,
                "e_corr": row.e_corr,
                "e_total": row.e_total,
                "e_corr_eV": row.e_corr * EH_EV,
                "e_total_eV": row.e_total * EH_EV,
                "l2": _expect(ops["L2"], ci_p, ncas, na, nb),
                "lz": _expect(ops["Lz"], ci_p, ncas, na, nb),
                "s2": _expect(ops["S2"], ci_p, ncas, na, nb),
                "classes": [{
                    "label": c.label, "note": c.note,
                    "norm": c.norm, "energy": c.energy,
                    "energy_meh": c.energy * 1e3,
                    "e_per_norm": c.e_per_norm,
                } for c in row.classes],
            })
        print(f"  [rows] {addr}: nroots {nroots}, "
              f"E(2) root0 {rows[0].e_corr * 1e3:+.4f} mEh "
              f"({rows[0].e_corr * EH_EV:+.4f} eV), "
              f"<L2> {entry['roots'][0]['l2']:.6f} "
              f"({time.time() - t0:.1f} s)", flush=True)
        entries.append(entry)
    return {
        "atom": atom,
        "provenance": {
            "engine": "nevpt2-sc (pyscf.mrpt.NEVPT, strongly contracted)",
            "pyscf_version": pyscf_version,
            "collected": datetime.now().isoformat(timespec="seconds"),
            "date": date.today().isoformat(),
            "declared_roots_policy": "every stored root is a row; no "
                                     "silent draws",
            "e_ref_partition": "kernel totals only (the readout's own "
                               "root totals); e_total = e_ref + e_corr",
            "wobble_note": "SC-NEVPT2 wobbles ~3e-7 Eh run-to-run "
                           "(thread-parallel reductions); the committed "
                           "payload is one draw, the store's chk is the "
                           "reproducibility unit",
            "stored_roots_bands": "energies 1e-11 absolute; CI "
                                  "amplitudes 1e-8 per-root; the "
                                  "degenerate block's joint projector; "
                                  "the re-solve angle (davidson noise) "
                                  "1e-3",
            "meters": "L2/Lz/S2 measured live through the CI-order "
                      "bridge on the stored draw",
        },
        "cells": entries,
    }


def main() -> None:
    atoms = sys.argv[1:]
    if not atoms:
        sys.exit("usage: run_pt2_ladder.py <atom> [<atom> ...]")
    for atom in atoms:
        print(f"===== pt2 ladder: {atom} =====", flush=True)
        payload = collect_atom(atom)
        out = Path("elements") / atom / "perturbation.json"
        out.write_text(json.dumps(payload, indent=1) + "\n")
        n_rows = sum(c["disposition"] == "rows" for c in payload["cells"])
        n_ref = len(payload["cells"]) - n_rows
        print(f"wrote {out}: {n_rows} rows / {n_ref} refused "
              f"of {len(payload['cells'])} declared cells", flush=True)


if __name__ == "__main__":
    main()
