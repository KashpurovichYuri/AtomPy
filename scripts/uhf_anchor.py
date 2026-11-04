"""Warm the UHF anchor record into an atom's cells.json: the gated
run_uhf (:attr:`ReferenceSpec.occ_e_ref_uhf` -- the census measured
it, :func:`atomterm.reference.uhf_basin_census`) plus the
contamination reading.  The ladder builder renders the ``UHF
(anchor)`` row from this record; re-running is idempotent (the
cached chk is trusted, the record is overwritten with the same
numbers).

Usage:
    python scripts/uhf_anchor.py <Atom>
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from atomterm.cases import CAMPAIGN  # noqa: E402
from atomterm.reference import run_uhf, uhf_contamination  # noqa: E402


def main(atom: str) -> None:
    spec = CAMPAIGN[atom]
    if spec.ref is None or spec.ref.occ_e_ref_uhf is None:
        sys.exit(f"{atom}: no UHF gate declared "
                 "(ReferenceSpec.occ_e_ref_uhf) -- the census runs "
                 "first, the declaration picks the basin")
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
