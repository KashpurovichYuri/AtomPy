"""The frozen-scan anatomy runner: one scan report per declared atom.

The scan data (families, angle grids, fixed extras) is DECLARED --
``CaseSpec.scan`` in ``cases.CAMPAIGN``; the machinery is
:mod:`atomterm.scans` (optimizer-free, deterministic given the ROHF
reference -- the campaign cache ``elements/<atom>/rungs/rohf.chk``).
Writes ``elements/<atom>/scan.md`` and prints it.

Usage:
    python scripts/radial_scan.py            # every atom declaring a scan
    python scripts/radial_scan.py C Sc       # the named atoms only
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402

from atomterm.cases import CAMPAIGN  # noqa: E402
from atomterm.reference import run_rohf  # noqa: E402
from atomterm.scans import member_tilt_map, radial_scan, scan_state, \
    scan_window  # noqa: E402
from atomterm.windows import select_shell_ao  # noqa: E402


def run_one(name: str, spec) -> None:
    """One atom's scan: the m-symmetric curve over the declared grid,
    then the member-tilt map at its optimum; rendered to stdout and
    ``elements/<atom>/scan.md``."""
    cfg = spec.atom_config()
    chk = Path("elements") / name / "rungs" / "rohf.chk"
    chk.parent.mkdir(parents=True, exist_ok=True)
    mol, mf = run_rohf(cfg, str(chk), ref=spec.ref)

    fam_a = select_shell_ao(mol, spec.scan.fam[0])
    fam_b = select_shell_ao(mol, spec.scan.fam[1])
    ncore = (mol.nelectron - spec.nelec_act) // 2
    na = (spec.nelec_act + spec.spin) // 2
    nb = (spec.nelec_act - spec.spin) // 2

    extra = None
    if spec.scan.extra_fixed:
        col = np.eye(mol.nao_nr())[:, select_shell_ao(mol,
                                                      spec.scan.extra_fixed)]
        extra = col.reshape(-1, 1)

    def state(thetas):
        c_act = scan_window(mol, mf, thetas, fam_a, fam_b)
        if extra is not None:
            c_act = np.hstack([c_act, extra])
        return scan_state(mol, mf, c_act, ncore, na, nb)

    lines = [f"# Frozen radial scan: {name} "
             f"({spec.scan.fam[0]} -> {spec.scan.fam[1]})",
             "",
             "m-symmetric curve (all members share the angle: EXACTLY",
             "closed windows, [H, Lz] = 0 throughout):", "",
             f"{'theta':>7} {'dE(uEh)':>11} {'<H1>+C(uEh)':>13} "
             f"{'<H2>(uEh)':>11} {'|[H,Lz]|':>9}"]
    e_zero = None
    best = None
    best_th = None
    for th in np.linspace(*spec.scan.grid):
        p = state([th] * len(fam_a))
        if e_zero is None:
            e_zero = p.e_active
        if best is None or p.e_active < best.e_active:
            best, best_th = p, th
        lines.append(f"{th:7.3f} {(p.e_active - e_zero) * 1e6:11.1f} "
                     f"{p.e_one * 1e6:13.1f} {p.e_two * 1e6:11.1f} "
                     f"{p.comm_h_lz:9.1e}")

    lines += ["", f"optimum: theta* = {best_th:.3f}, "
              f"dE = {(best.e_active - e_zero) * 1e6:+.1f} uEh "
              f"(compare the free/radial CASSCF gains, AGENTS.md fact 16), "
              f"|[H,Lz]| = {best.comm_h_lz:.1e}", "",
              f"member-tilt map at theta* = {best_th:.3f} "
              f"(delta = 0.10 rad):", "",
              f"{'tilted':>12} {'dE(uEh)':>11} {'|[H,Lz]|':>9}"]
    tilts = member_tilt_map(mol, mf, ncore, na, nb, fam_a, fam_b,
                            best_th, 0.10, extra_cols=extra)
    for key in sorted(tilts, key=lambda k: len(k)):
        p = tilts[key]
        lines.append(f"{key:>12} {(p.e_active - best.e_active) * 1e6:11.2f} "
                     f"{p.comm_h_lz:9.2e}")

    text = "\n".join(lines)
    print(text)
    out = Path("elements") / name / "scan.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text + "\n", encoding="utf-8")
    print(f"\nwrote {out}")


def main() -> None:
    """Every named atom's declared scan (default: all declared scans)."""
    atoms = sys.argv[1:] or [a for a, s in CAMPAIGN.items() if s.scan]
    unknown = [a for a in atoms if a not in CAMPAIGN]
    if unknown:
        sys.exit(f"not declared: {', '.join(unknown)}; "
                 f"declared: {', '.join(CAMPAIGN)}")
    for name in atoms:
        if CAMPAIGN[name].scan is None:
            print(f"{name}: no declared scan -- skipped")
            continue
        run_one(name, CAMPAIGN[name])


if __name__ == "__main__":
    main()
