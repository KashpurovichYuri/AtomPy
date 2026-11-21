"""The U14e MOM census: the members under the ranking-free follow.

Step 0's attribution (elements/U14e/logs/spectrum.log): the weight
hook's failure on this atom is a VOCABULARY defect -- the premium
sees only character, and ~12 pure d/f valence virtuals plus the
tight-primitive tail all read premium ~ 2.000 (weights exactly 1.000,
no metric noise).  The certified cure class is MOM: the occupations
FOLLOW a seed by max S-overlap instead of ranking
(reference._mom_occ -- the twins' repair; on a metric with NO
near-null tail the follow's guarantee is cleaner still).

This census arms _mom_occ with the head's OWN columns re-assigned to
each configuration (no foreign source basis exists for U -- this is
the self-seeded form the seeded mode cannot express).  The seed
measurement decides the seam: if the members land here, the
explicit-seed master joins run_rohf's vocabulary and the series
revives; if not, the refusal stands with the landscape attribution.

Usage:
    OMP_NUM_THREADS=16 python scripts/probe_u14e_mom.py \
        > elements/U14e/logs/mom_census.log
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from atomterm.cases import campaign_spec  # noqa: E402
from atomterm.reference import (_accelerated, _mom_occ,  # noqa: E402
                                occupation_postcheck, run_rohf)
from atomterm.windows import pattern_projector  # noqa: E402

MEMBERS = (
    ("f4", {"l3": 4}),
    ("f2d2", {"l3": 2, "l2": 2}),
)
ACCELERATORS = ("diis", "diis+newton", "newton")
LS = ("l0", "l1", "l2", "l3")


def fingerprint(mol, mf, projs):
    opens = np.flatnonzero(np.asarray(mf.mo_occ) == 1.0)
    c = mf.mo_coeff[:, opens]
    return {p: float(np.einsum("ai,ab,bi->", c, P, c).real)
            for p, P in projs.items()}


def fmt(fp):
    return "  ".join(f"{k}:{v:.3f}" for k, v in fp.items())


def seed_columns(mf, want):
    """The head's columns re-assigned to a configuration: the closed
    block untouched; the open block re-picked by character (the
    reference's open columns of that l deepest-first, then its best
    virtual columns of the same l -- highest weight, lowest epsilon)
    until the count is met.  Returns (c_seed, occ_seed)."""
    occ = np.asarray(mf.mo_occ)
    c, eps = mf.mo_coeff, np.asarray(mf.mo_energy)
    closed = np.flatnonzero(occ == 2.0)
    opens = np.flatnonzero(occ == 1.0)
    virt = np.flatnonzero(occ == 0.0)
    projs = {p: pattern_projector(mf.mol, (p,)) for p in LS}
    pick = []
    for p, n_want in want.items():
        P = projs[p]

        def w(cols):
            return np.einsum("ai,ab,bi->i", c[:, cols].conj(), P,
                             c[:, cols]).real

        have = sorted((i for i in opens if w([i])[0] >= 0.5),
                      key=lambda i: eps[i])[:n_want]
        pick += have
        short = n_want - len(have)
        if short > 0:
            # EPSILON-FIRST within the majority floor (w >= 0.5), never
            # weight-ranked: at full precision the tight-primitive
            # virtuals' weights edge out the valence members (the f4
            # seed's 4th pick read eps = 22.4 under a (-w, eps) sort --
            # the diagnosis's own defect, the character ranking cannot
            # tell a 5f member from a 14f primitive).  The shell is
            # selected by energy among the columns that carry the
            # wanted character
            cand = sorted((i for i in virt if w([i])[0] >= 0.5),
                          key=lambda i: eps[i])[:short]
            if len(cand) < short:
                raise RuntimeError(
                    f"the seed wants {n_want} {p}-character columns, "
                    "the reference does not carry them")
            pick += cand
    cols = np.concatenate([closed, np.asarray(pick, dtype=int)])
    occ_seed = np.concatenate([np.full(len(closed), 2.0),
                               np.ones(len(pick))])
    return c[:, cols], occ_seed


def main() -> None:
    spec = campaign_spec("U14e")
    cfg = spec.atom_config()
    scratch = Path("elements/U14e/rungs/mom_tmp")
    scratch.mkdir(parents=True, exist_ok=True)
    for f in scratch.glob("*.chk"):        # a killed prior run's
        f.unlink()                         # scratch must not cache
    mol, mf = run_rohf(cfg, "elements/U14e/rungs/rohf.chk", ref=spec.ref)
    projs = {p: pattern_projector(mol, (p,)) for p in LS}
    print(f"head (trusted): E = {mf.e_tot:.8f}  open-block "
          f"{fmt(fingerprint(mol, mf, projs))}")

    for tag, want in MEMBERS:
        c_seed, occ_seed = seed_columns(mf, want)
        dm0 = (c_seed * occ_seed) @ c_seed.conj().T
        dm0 = dm0.real.astype(np.float64)
        print(f"\n== member {tag}  seed want {want}  "
              f"(MOM follows the re-assigned head columns)")
        for acc in ACCELERATORS:
            m = __import__("pyscf").scf.ROHF(mol)
            m.get_occ = _mom_occ(mol, c_seed, occ_seed)
            m.conv_tol = 1e-12
            from atomterm.reference import ReferenceSpec
            mf_x = _accelerated(m, ReferenceSpec(conv_accelerator=acc),
                                dm0=dm0)
            fp = fingerprint(mol, mf_x, projs)
            print(f"  MOM {acc:>12}: E = {mf_x.e_tot:.8f} "
                  f"conv={bool(mf_x.converged)}  | {fmt(fp)}", flush=True)
            if mf_x.converged:
                # the certified repeat arm: the same solution under a
                # fresh MOM start (stability, not a one-off draw)
                m2 = __import__("pyscf").scf.ROHF(mol)
                m2.get_occ = _mom_occ(mol, mf_x.mo_coeff,
                                      np.asarray(mf_x.mo_occ))
                m2.conv_tol = 1e-12
                m2.mo_coeff = mf_x.mo_coeff
                m2.mo_energy = mf_x.mo_energy
                m2.mo_occ = mf_x.mo_occ
                mf2 = _accelerated(m2, ReferenceSpec(
                    conv_accelerator="diis"))
                print(f"  {'  repeat (diis)':>18}: E = {mf2.e_tot:.8f} "
                      f"conv={bool(mf2.converged)}  | "
                      f"{fmt(fingerprint(mol, mf2, projs))}", flush=True)

    for f in scratch.glob("*.chk"):
        f.unlink()
    scratch.rmdir()
    print("\ndone -- a landed member is a converged MOM arm whose "
          "fingerprint carries its configuration")


if __name__ == "__main__":
    main()
