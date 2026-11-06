"""The measured shell composition of a converged reference -- the
certificate behind a declared shell layer (``CaseSpec.shell_bands``).

On an uncontracted basis the AO labels are radial counters; which
bands realize a chemical shell is a MEASURED fact (this script) that
the declaration then states (the layer) and the folds/windows consume
(:func:`atomterm.windows.shell_realization`).

TWO METERS, one per role:
- the DISPLAY meter: per-band Mulliken populations of the canonical
  MOs, ``w_b(c) = sum_{i in b} c_i (S c)_i`` -- partitioning unity per
  MO exactly (each AO row belongs to one band).  On near-dependent
  primitive bands these run WILD (large positives and negatives
  cancelling); they are the eyeball table, never the certificate.
- the CERTIFICATE meter: the S-projector capture ``|C_span^dag S c|^2``
  (bounded [0, 1] whatever the band overlaps) -- the same weight
  :func:`atomterm.core.confined_core_orbitals` ranks and records.  The
  suggested realization per chemical shell is the GREEDY FORWARD
  band set (each step adding the band that most lifts the group's
  WORST-MO capture) until every MO of the shell is captured at
  --target.

The chemical grouping is node-counting (the k-th l-type occupied MO by
energy is the k-th n of that l).  The SCF here is a PLAIN substrate
(any converged filling): the deep closed MOs are basin-independent,
and only their band capture matters.  The DECLARATION is the caller's
decision; this log is its evidence.

Usage:
    python scripts/band_composition.py Sc data/sc_ccpvdz_unc.txt --spin 1 \
        > elements/Sc-unc/logs/band_composition.log
"""

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pyscf import gto, scf  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from atomterm.windows import lowdin, select_shell_ao  # noqa: E402

_FIRST_N = {"s": 1, "p": 2, "d": 3, "f": 4, "g": 5}
_ANG = "spdfgh"


def band_labels(mol) -> list[str]:
    """The AO band label per AO row (a band = one label token with the
    Cartesian member suffixes stripped -- a band is the full m-set,
    the unit the shell layer declares)."""
    import re

    return [re.sub(r"[xyz0-9^+-]*$", "", line.split()[-1])
            for line in mol.ao_labels()]


def band_populations(mol, mo_coeff) -> dict[str, np.ndarray]:
    """Per-band Mulliken populations of every MO column (the DISPLAY
    meter; see the module docstring for the certificate meter)."""
    s = mol.intor("int1e_ovlp")
    sc = s @ mo_coeff                            # (nao, nmo): (S c)[j, m]
    labels = band_labels(mol)
    out: dict[str, np.ndarray] = {}
    for b in dict.fromkeys(labels):             # first-appearance order
        rows = [i for i, lbl in enumerate(labels) if lbl == b]
        out[b] = np.einsum("im,im->m", mo_coeff[rows], sc[rows]).real
    tot = np.sum([w for w in out.values()], axis=0)
    assert np.abs(tot - 1.0).max() < 1e-9, "band populations must sum to 1"
    return out


def capture(s, span_rows: list[int], cols: np.ndarray) -> np.ndarray:
    """The CERTIFICATE meter: ``|C^dag S c|^2`` per column of ``cols``
    on the Loewdin span of the AO rows (bounded [0, 1])."""
    c_span = lowdin(np.eye(s.shape[0])[:, span_rows], s)
    ov = c_span.conj().T @ s @ cols
    return np.einsum("mi,mi->i", ov, ov.conj()).real


def greedy_realization(mol, s, mo_coeff, mos: list[int], shell: str,
                       target: float) -> tuple[list[str], float]:
    """The greedy forward band set realizing one chemical shell: each
    step adds the band (of the shell's angular momentum) that most
    lifts the group's WORST-MO projector capture, until every MO of
    the group is captured at ``target`` (or the band supply runs out).

    Returns ``(bands, worst_capture)`` in selection order.
    """
    labels = band_labels(mol)
    bands = [b for b in dict.fromkeys(labels) if b[-1] == shell[-1]]
    rows_of = {b: select_shell_ao(mol, b) for b in bands}
    cols = np.asarray(mo_coeff)[:, mos]
    picked: list[str] = []
    worst = 0.0
    while len(picked) < len(bands):
        best_b, best_w = None, -1.0
        for b in bands:
            if b in picked:
                continue
            rows = sorted(sum((rows_of[p] for p in picked + [b]), []))
            try:
                w = capture(s, rows, cols).min()
            except ValueError:                   # near-dependent set
                continue
            if w > best_w:
                best_b, best_w = b, w
        if best_b is None:
            break
        picked.append(best_b)
        worst = best_w
        if worst >= target:
            break
    return picked, worst


def main(argv) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("symbol")
    p.add_argument("basis", help="basis label or file path")
    p.add_argument("--spin", type=int, default=0)
    p.add_argument("--ecp", default=None)
    p.add_argument("--target", type=float, default=0.995,
                   help="the worst-MO projector capture a suggested "
                        "realization must reach")
    p.add_argument("--floor", type=float, default=5e-3,
                   help="per-MO Mulliken weight above which a band "
                        "prints in the display table")
    args = p.parse_args(argv[1:])

    basis = str(Path(args.basis)) if Path(args.basis).exists() else args.basis
    mol = gto.M(atom=f"{args.symbol} 0 0 0", basis=basis, spin=args.spin,
                ecp=args.ecp, verbose=0)
    mf = scf.ROHF(mol)
    mf.conv_tol = 1e-11
    mf.kernel()
    if not mf.converged:
        print(f"# plain DIIS did not converge (E = {mf.e_tot:.6f}) "
              "-- retrying Newton", file=sys.stderr)
        mf = mf.newton()
        mf.kernel()
    if not mf.converged:
        sys.exit("the measurement SCF did not converge")
    print(f"# {args.symbol} plain ROHF E = {mf.e_tot:.8f} Eh "
          f"(measurement substrate; the closed MOs are basin-independent)")

    mo_coeff = np.asarray(mf.mo_coeff)
    mo_occ = np.asarray(mf.mo_occ)
    mo_energy = np.asarray(mf.mo_energy)
    pops = band_populations(mol, mo_coeff)
    s = mol.intor("int1e_ovlp")

    # the chemical grouping: within l, a closed shell of an atom is a
    # COMPLETE m-set (2l+1 MOs) -- chunk the closed MOs by
    # multiplicity in energy order (the node-counted n-sequence);
    # C1 degeneracy gaps are numerically ragged, so raw-energy
    # clustering is the fragile way.  The frontier's singly-occupied
    # MOs form the open shell (its partial m-set -- 3d1)
    dom_l = [max(_ANG, key=lambda l: sum(
        w[i] for b, w in pops.items() if b[-1] == l))
        for i in range(len(mo_occ))]
    groups: dict[str, list[int]] = {}
    for l in _ANG:
        m_set = 2 * _ANG.index(l) + 1
        closed = sorted(i for i, (o, d) in enumerate(zip(mo_occ, dom_l))
                        if d == l and o > 1.5)
        open_ = sorted(i for i, (o, d) in enumerate(zip(mo_occ, dom_l))
                       if d == l and 0.5 < o <= 1.5)
        if len(closed) % m_set:
            print(f"# WARNING: {len(closed)} closed {l}-MOs do not "
                  "chunk into whole m-sets -- the dominant-l reading "
                  "misfired somewhere", file=sys.stderr)
        chunks = [closed[i:i + m_set] for i in range(0, len(closed), m_set)]
        if open_:
            chunks.append(open_)
        for k, cl in enumerate(chunks):
            groups.setdefault(f"{_FIRST_N[l] + k}{l}", []).extend(cl)

    print(f"\n# the display table (Mulliken populations -- wild on "
          "near-dependent bands; the certificate is the projector "
          f"capture below; floor {args.floor:.1%})")
    for shell, mos in groups.items():
        wb = {b: sum(w[i] for i in mos) / len(mos) for b, w in pops.items()}
        top = [b for b in sorted(wb, key=lambda b: -wb[b])
               if wb[b] >= args.floor]
        print(f"{shell:>5} occ={mo_occ[mos[0]]:.0f} MOs={mos}  "
              + " ".join(f"{b}:{wb[b]:.2f}" for b in top))

    print(f"\n# the certificate: the greedy realization (worst-MO "
          f"projector capture >= {args.target:.1%}; the DECLARATION "
          "takes these lists)")
    for shell, mos in groups.items():
        bands, worst = greedy_realization(mol, s, mo_coeff, mos, shell,
                                          args.target)
        note = "" if worst >= args.target else (
            f"  !! supply exhausted at {worst:.3f}")
        print(f"{shell:>5} -> ({', '.join(bands)})  capture {worst:.4f}"
              f"{note}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
