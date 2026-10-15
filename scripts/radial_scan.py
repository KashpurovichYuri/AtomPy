"""Frozen radial scans for one atom (ticket 16, Phase 1): the
m-symmetric correlation curve at exact symmetry + the member-tilt map.

Usage:
    python scripts/radial_scan.py [C|Si|Sc] [workdir]

Writes results/<label>_scan.md and prints it. No optimizer anywhere --
every window is an explicit cos/sin construction (see atomterm.scans).
"""

import sys
from pathlib import Path

import numpy as np

from atomterm.drivers import run_rohf, select_shell_ao
from atomterm.scans import member_tilt_map, radial_scan, scan_state, scan_window
from atomterm.system import AtomConfig

SPECS = {
    "C": dict(cfg=AtomConfig(), fam=("2p", "3p"), nelec=2, spin=2,
              grid=np.linspace(0.0, 0.8, 17)),
    "Si": dict(cfg=AtomConfig("Si", "ccpvdz", spin=2), fam=("3p", "4p"),
               nelec=2, spin=2, grid=np.linspace(0.0, 1.1, 23)),
    # d-only tilts (4s column appended unscanned: the 4s->5s direction
    # is harmful and would mask the d structure)
    "Sc": dict(cfg=AtomConfig("Sc", "ccpvdz", spin=1), fam=("3d", "4d"),
               nelec=3, spin=1, grid=np.linspace(0.0, 0.4, 9),
               extra_fixed="4s"),
}


def main():
    sym = sys.argv[1] if len(sys.argv) > 1 else "C"
    workdir = Path(sys.argv[2] if len(sys.argv) > 2 else "results")
    workdir.mkdir(parents=True, exist_ok=True)
    spec = SPECS[sym]
    cfg = spec["cfg"]
    mol, mf = run_rohf(cfg, str(workdir / f"{cfg.label}_rohf.chk"))

    fam_a = select_shell_ao(mol, spec["fam"][0])
    fam_b = select_shell_ao(mol, spec["fam"][1])
    ncore = (mol.nelectron - spec["nelec"]) // 2
    na = (spec["nelec"] + spec["spin"]) // 2
    nb = (spec["nelec"] - spec["spin"]) // 2

    extra = None
    if spec.get("extra_fixed"):
        col = np.eye(mol.nao_nr())[:, select_shell_ao(mol, spec["extra_fixed"])]
        extra = col.reshape(-1, 1)

    def state(thetas):
        c_act = scan_window(mol, mf, thetas, fam_a, fam_b)
        if extra is not None:
            c_act = np.hstack([c_act, extra])
        return scan_state(mol, mf, c_act, ncore, na, nb)

    lines = [f"# Frozen radial scan: {sym} ({spec['fam'][0]} -> {spec['fam'][1]})",
             "",
             "m-symmetric curve (all members share the angle: EXACTLY",
             "closed windows, [H, Lz] = 0 throughout):", "",
             f"{'theta':>7} {'dE(uEh)':>11} {'<H1>+C(uEh)':>13} {'<H2>(uEh)':>11} {'|[H,Lz]|':>9}"]
    e_zero = None
    best = None
    best_th = None
    for th in spec["grid"]:
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
             f"member-tilt map at theta* = {best_th:.3f} (delta = 0.10 rad):", "",
             f"{'tilted':>12} {'dE(uEh)':>11} {'|[H,Lz]|':>9}"]
    tilts = member_tilt_map(mol, mf, ncore, na, nb, fam_a, fam_b,
                            best_th, 0.10, extra_cols=extra)
    for key in sorted(tilts, key=lambda k: len(k)):
        p = tilts[key]
        lines.append(f"{key:>12} {(p.e_active - best.e_active) * 1e6:11.2f} "
                     f"{p.comm_h_lz:9.2e}")

    text = "\n".join(lines)
    print(text)
    (workdir / f"{cfg.label}_scan.md").write_text(text + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
