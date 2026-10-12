#!/usr/bin/env python
"""Run the carbon term-symbol pipeline: ROHF -> CASCI -> CASSCF -> SA-CASSCF.

Usage:
    python scripts/run_carbon.py [workdir]
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from atomterm.pipeline import run_carbon  # noqa: E402

if __name__ == "__main__":
    workdir = sys.argv[1] if len(sys.argv) > 1 else "results"
    run_carbon(workdir=workdir)
