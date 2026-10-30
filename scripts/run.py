#!/usr/bin/env python
"""Run the declared campaign: the term tables per atom.

Each named atom executes its :class:`atomterm.cases.CaseSpec` (ROHF ->
CASCI -> CASSCF -> SA-CASSCF, checkpointed under ``elements/<Atom>/``);
an atom is runnable iff declared -- new atoms enter through
``cases.CAMPAIGN``, never through this script. The methods grid and
its tables are the cells deliverable (``scripts/regen_cells.py`` +
the two renderers), not this one.

Usage:
    python scripts/run.py            # every declared atom, CAMPAIGN order
    python scripts/run.py Ti C       # the named atoms only
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from atomterm.cases import CAMPAIGN, run_case  # noqa: E402


def main():
    atoms = sys.argv[1:] or list(CAMPAIGN)
    unknown = [a for a in atoms if a not in CAMPAIGN]
    if unknown:
        sys.exit(f"not declared: {', '.join(unknown)}; "
                 f"declared: {', '.join(CAMPAIGN)}")
    for atom in atoms:
        run_case(CAMPAIGN[atom])


if __name__ == "__main__":
    main()
