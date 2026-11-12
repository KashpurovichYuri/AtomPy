"""STAGE 1a (ticket 37): the overlap-spectrum measurement on the
uncontracted twins -- the diagnosis the conditioning design rests on.

Per l-block eigenvalues of S (a single atom, spherical basis: S is
block-diagonal in l, so the per-l block is exact).  The GO/NO-GO fact:
does the spectrum show a GAP (meaningful ~1e-2 vs null ~1e-7), or a
continuous decay?

Sandbox protocol (ticket 37): in-repo, imported by nothing, deleted at
integration.  Run: OMP_NUM_THREADS=64 PYTHONPATH=src python
sandbox/spectrum.py > elements/<twin>/logs/spectrum.log (per atom, or
both to stdout).
"""
import numpy as np

from atomterm.cases import campaign_spec
from atomterm.system import build_mol

# the thresholds the design speaks in (the gap vocabulary: meaningful
# vs null) -- print the counts around each
TOL = (1e-1, 1e-2, 1e-4, 1e-6, 1e-8, 1e-10)


def l_blocks(mol):
    """AO index lists per angular momentum (a single atom: exact)."""
    out = {}
    ao_loc = mol.ao_loc_nr()
    for i in range(mol.nbas):
        l = mol.bas_angular(i)
        p0, p1 = int(ao_loc[i]), int(ao_loc[i + 1])
        out.setdefault(l, []).extend(range(p0, p1))
    return out


def spectrum(key: str) -> str:
    cfg = campaign_spec(key).atom_config()
    mol = build_mol(cfg)
    s = mol.intor("int1e_ovlp")
    lines = [f"{key}  ({cfg.basis if isinstance(cfg.basis, str) else key}"
             f" basis; nao = {mol.nao})"]
    for l, ix in sorted(l_blocks(mol).items()):
        ev = np.linalg.eigvalsh(s[np.ix_(ix, ix)])[::-1]
        n = len(ev)
        counts = " ".join(f">{t:.0e}:{sum(ev > t)}" for t in TOL)
        lines.append(f"  l={l}  n={n}  {counts}")
        # the gap region: every eigenvalue below 1e-1, descending
        tail = [e for e in ev if e < 1e-1]
        lines.append("    tail: " + " ".join(f"{e:.2e}" for e in tail))
        # the maximal ratio jump in the tail = the gap candidate
        if len(tail) > 1:
            ratios = [tail[i] / tail[i + 1] for i in range(len(tail) - 1)]
            j = int(np.argmax(ratios))
            lines.append(f"    max gap: {tail[j]:.2e} -> {tail[j+1]:.2e}"
                         f" (x{ratios[j]:.0f}, {j+1} above)")
    return "\n".join(lines)


if __name__ == "__main__":
    for key in ("Sc-unc", "Ti-unc"):
        text = spectrum(key)
        print(text)
        import pathlib
        pathlib.Path(f"elements/{key}/logs").mkdir(parents=True,
                                                   exist_ok=True)
        pathlib.Path(f"elements/{key}/logs/spectrum.log").write_text(
            text + "\n")
