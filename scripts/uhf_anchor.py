"""Warm the UHF anchor record into an atom's cells.json: run_uhf
(:attr:`ReferenceSpec.occ_e_ref_uhf` gates it where the census shows
multiple basins -- the census measures, the declaration picks; a
single-basin atom runs ungated, the census verdict being its guard)
plus the contamination reading.  The methods-table builder renders the
``UHF (anchor)`` row from this record; re-running is idempotent (the
cached chk is trusted, the record is overwritten with the same
numbers).

Usage:
    python scripts/uhf_anchor.py <Atom>
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from atomterm.cases import campaign_spec  # noqa: E402
from atomterm.reference import run_uhf, uhf_contamination  # noqa: E402


def main(atom: str) -> None:
    spec = campaign_spec(atom)
    store = Path("elements") / atom / "cells.json"
    data = json.loads(store.read_text())
    chk = str(Path("elements") / atom / "rungs" / "uhf.chk")
    _, mf = run_uhf(spec.atom_config(), chk, ref=spec.ref)
    s2, sigma = uhf_contamination(mf)
    data["uhf"] = {"e_tot": float(mf.e_tot), "s2": float(s2),
                   "sigma_min": float(sigma.min()),
                   "converged": bool(mf.converged)}
    store.write_text(json.dumps(data, indent=2))
    print(f"{atom}: UHF E = {mf.e_tot:.8f} Eh, <S^2> = {s2:.6f}, "
          f"sigma_min = {float(sigma.min()):.6f} -> {store}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "Sc")
