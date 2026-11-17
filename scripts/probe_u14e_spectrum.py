"""The U14e basis's conditioning diagnosis (ticket 39's follow-up wave,
Step 0): the per-l overlap spectra, the amplification the member hunts
see, and the distillation chains' margins.

The instrument is the twins' Stage-1a port (git history,
sandbox/spectrum.py at 2283f3d) with three additions the U14e question
needs:

- kappa per l-block: the twins' amplification measured kappa ~ 2000
  (moderate conditioning, NOT near-linear dependence) -- and their
  cut was NO-GO even so (no gap).  The cut decision here needs the
  same separation: a kappa statement is about amplification, a gap
  statement is about cutting.
- the tail's location: which band primitives dominate the smallest
  eigenvectors (diffuse vs tight end) -- where a cut, if ever
  warranted, would live.
- the hook-visible garbage: the trusted head's virtual columns that
  carry l-span weights ~1 at positive epsilon -- the amplified set
  the weight hook's premium ranks against (the member hunts' failure
  signature, quantified on THIS basis).

Usage:
    OMP_NUM_THREADS=8 python scripts/probe_u14e_spectrum.py \
        > elements/U14e/logs/spectrum.log
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from atomterm.cases import campaign_spec  # noqa: E402
from atomterm.reference import run_rohf  # noqa: E402
from atomterm.windows import member_chains, pattern_projector  # noqa: E402

# the twins' gap vocabulary (Stage 1a): meaningful ~1e-2 vs null ~1e-7
TOL = (1e-1, 1e-2, 1e-4, 1e-6, 1e-8, 1e-10)
LETTERS = "spdfgh"


def l_blocks(mol):
    """AO index lists per angular momentum (a single atom: exact)."""
    out = {}
    ao_loc = mol.ao_loc_nr()
    for i in range(mol.nbas):
        l = mol.bas_angular(i)
        p0, p1 = int(ao_loc[i]), int(ao_loc[i + 1])
        out.setdefault(l, []).extend(range(p0, p1))
    return out


def band_labels(mol, ix):
    """The band label per AO row of one l-block (shell counters -- the
    honest tag for locating the tail)."""
    labels = [line.split()[-1] for line in mol.ao_labels()]
    return [labels[i] for i in ix]


def spectrum_block(mol, s, l, ix):
    """One l-block's spectrum, kappa, tail location, and gap read."""
    ev, vec = np.linalg.eigh(s[np.ix_(ix, ix)])
    order = np.argsort(ev)[::-1]
    ev, vec = ev[order], vec[:, order]
    n = len(ev)
    kappa = ev[0] / ev[-1] if ev[-1] > 0 else float("inf")
    lines = [f"  l={l} ({LETTERS[l]})  n={n}  "
             + " ".join(f">{t:.0e}:{sum(ev > t)}" for t in TOL)
             + f"  kappa={kappa:.2e}  min={ev[-1]:.2e}  max={ev[0]:.2e}"]
    tail_mask = ev < 1e-1
    tail = ev[tail_mask]
    if len(tail):
        lines.append("    tail: " + " ".join(f"{e:.2e}" for e in tail))
        if len(tail) > 1:
            ratios = tail[:-1] / tail[1:]
            j = int(np.argmax(ratios))
            lines.append(f"    max gap: {tail[j]:.2e} -> {tail[j+1]:.2e}"
                         f" (x{ratios[j]:.0f}, {j+1} above)")
        # the location: the three smallest eigenvectors' weights per band
        bands = band_labels(mol, ix)
        loc: dict[str, float] = {}
        for v in vec[:, -3:].T:
            for b, w in zip(bands, np.abs(v)):
                loc[b] = loc.get(b, 0.0) + float(w) ** 2
        top = sorted(loc.items(), key=lambda kv: -kv[1])[:4]
        lines.append("    smallest-3 eigvec weight: "
                     + ", ".join(f"{b}:{w:.2f}" for b, w in top))
    else:
        lines.append("    tail: (none below 1e-1)")
    return lines


def main() -> None:
    spec = campaign_spec("U14e")
    cfg = spec.atom_config()
    mol, mf = run_rohf(cfg, "elements/U14e/rungs/rohf.chk", ref=spec.ref)
    s = mol.intor("int1e_ovlp")
    nao = mol.nao_nr()
    print(f"U14e  (nao = {nao})")
    blocks = l_blocks(mol)
    for l, ix in sorted(blocks.items()):
        for line in spectrum_block(mol, s, l, ix):
            print(line)

    # the distillation chains' margins (what metric_radials' Cholesky
    # inverts): the prior families' s_r blocks
    print("\n== prior-chain conditioning (the distillations' own s_r)")
    for shell, bands in (("5f", ("9f", "10f", "11f")),
                         ("6d", ("12d", "13d", "14d"))):
        chain, s_r = member_chains(mol, shell, bands, s)
        ev = np.linalg.eigvalsh(s_r)[::-1]
        print(f"  {shell} {'+'.join(bands)}: s_r {len(ev)}x{len(ev)} "
              f"eigenvalues {[f'{e:.4f}' for e in ev]}  "
              f"kappa={ev[0] / ev[-1]:.2f}")

    # the hook-visible garbage: the trusted head's virtual columns
    # carrying l-span weights ~1 at positive epsilon (the amplified
    # set the premium ranks against; the twins' diag signature)
    print("\n== the head's virtual columns (the hook's ranking field)")
    occ = np.asarray(mf.mo_occ)
    eps = np.asarray(mf.mo_energy)
    c = mf.mo_coeff
    projs = {p: pattern_projector(mol, (p,)) for p in ("l2", "l3")}
    w = {p: np.einsum("ai,ab,bi->i", c.conj(), P, c).real
         for p, P in projs.items()}
    wdwf = w["l2"] + w["l3"]
    virt = np.flatnonzero(occ == 0.0)
    sig = virt[wdwf[virt] > 0.5]
    print(f"  virtual columns with d+f weight > 0.5: {len(sig)}; "
          f"epsilon range {eps[sig].min():.3f} .. {eps[sig].max():.2f}")
    hi = sig[eps[sig] > 0.1]
    print(f"  of which at epsilon > 0.1 (the premium's garbage field): "
          f"{len(hi)}")
    if len(hi):
        order = hi[np.argsort(eps[hi])]
        rows = [f"    col {i}: eps={eps[i]:9.2f}  w_d={w['l2'][i]:.3f} "
                f"w_f={w['l3'][i]:.3f}" for i in order[:12]]
        print("\n".join(rows) + ("\n    ..." if len(order) > 12 else ""))
    occ_i = np.flatnonzero(occ > 0)
    print(f"  occupied columns for contrast (eps, w_d, w_f):")
    for i in occ_i:
        print(f"    col {i}: eps={eps[i]:9.3f}  w_d={w['l2'][i]:.3f} "
              f"w_f={w['l3'][i]:.3f}")


if __name__ == "__main__":
    main()
