"""The U14e series census: the f4 / f2d2 members' basins, measured.

The declaration's discovery half (discovery-then-declare): the
configuration series' non-head members do not exist yet -- their basin
energies (``SeriesSpec.e_ref``) and their accelerator recipes are what
this probe measures.  The arms see exactly what the production seam
can reach: the plain census cycle (accelerator x the gate's own
guesses, the census pattern -- hook attached, ``_accelerated``) under
the member's weight hook, plus ``run_rohf``'s caller-supplied ``dm0=``
form for a configuration the hook cannot PREFER -- the f2d2 member
shares the head's hook (same spans, same counts, same spin), so the
f3d1/f2d2 split is the SCF landscape's own and the seed is the honest
way to ask for the other basin.  Every arm reports its configuration
fingerprint (the open block's per-l weight sums): the character
post-check alone cannot discriminate f3d1 from f2d2 -- the fingerprint
and the e_ref gate own that discrimination.

No checkpoint survives the probe (scratch under the atom's rungs,
removed at exit); the DECLARED basin's certification runs -- gated,
to the trusted series paths -- are the caller's next step, in the
product format ``run_series`` loads trusted.

Usage:
    OMP_NUM_THREADS=16 python scripts/probe_u14e_series.py \
        > elements/U14e/logs/series_discovery.log
"""

import sys
from pathlib import Path

import numpy as np
from pyscf import scf

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from atomterm.cases import campaign_spec  # noqa: E402
from atomterm.reference import (ROHF_CENSUS_GUESSES, ReferenceSpec,  # noqa: E402
                                _accelerated, _forced_rohf_occ,
                                occupation_postcheck, run_rohf)
from atomterm.system import build_mol  # noqa: E402
from atomterm.windows import pattern_projector  # noqa: E402

# the declared members: the closed block stays the case reference's
# (the 6s6p7s pairs -- 5 pairs by the electron counts); the open spans
# name the configuration, ALL members at the head's spin (Hund: the
# f4 and f2d2 manifolds open at S = 2, the head f3d1 is spin 4 already)
MEMBERS = (
    ("f4", ("l0", "l1"), ("l3",)),
    ("f2d2", ("l0", "l1"), ("l2", "l3")),
)
ACCELERATORS = ("newton", "diis+newton", "diis")

# the fingerprint's vocabulary: per-l weight sums over the open block
LS = ("l0", "l1", "l2", "l3")


def lprojs(mol):
    return {p: pattern_projector(mol, (p,)) for p in LS}


def fingerprint(mol, mf, projs):
    """The open block's per-l weight sums (the configuration's own
    signature: f3d1 reads f~3 d~1, f2d2 f~2 d~2, f4 f~4)."""
    opens = np.flatnonzero(np.asarray(mf.mo_occ) == 1.0)
    c = mf.mo_coeff[:, opens]
    return {p: float(np.einsum("ai,ab,bi->", c, P, c).real)
            for p, P in projs.items()}


def fmt(fp):
    return "  ".join(f"{k}:{v:.3f}" for k, v in fp.items())


def report(label, mol, m, closed, opens, projs):
    print(f"  {label:>20}: E = {m.e_tot:.8f} conv={bool(m.converged)} "
          f"occ_ok={occupation_postcheck(mol, m, closed, opens)}"
          f"  | {fmt(fingerprint(mol, m, projs))}", flush=True)


def plain_arm(mol, closed, opens, acc, guess):
    """One census arm: the weight hook attached, the accelerator
    policy applied (the production path's own construction)."""
    m = scf.ROHF(mol)
    m.get_occ = _forced_rohf_occ(mol, list(closed), list(opens))
    m.conv_tol = 1e-12
    if guess is not None:
        m.init_guess = guess
    return _accelerated(m, ReferenceSpec(occ_closed=closed,
                                         occ_open=opens,
                                         conv_accelerator=acc))


def seed_density(mf, want):
    """The caller-supplied seed (``run_rohf(dm0=)``'s form): the
    trusted reference's own MOs re-assigned to a configuration -- the
    closed block untouched, the open block re-picked by CHARACTER
    (``want`` = per-l electron counts, e.g. {"l3": 2, "l2": 2}): the
    reference's open columns of that l (deepest first), then its best
    virtual columns of the same l (highest weight, lowest epsilon)
    until the count is met."""
    occ = np.asarray(mf.mo_occ)
    c, eps = mf.mo_coeff, np.asarray(mf.mo_energy)
    closed = np.flatnonzero(occ == 2.0)
    opens = np.flatnonzero(occ == 1.0)
    virt = np.flatnonzero(occ == 0.0)
    projs = lprojs(mf.mol)
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
            cand = sorted(virt, key=lambda i: (-w([i])[0], eps[i]))
            cand = [i for i in cand if w([i])[0] >= 0.5][:short]
            if len(cand) < short:
                raise RuntimeError(
                    f"the seed wants {n_want} {p}-character columns, "
                    "the reference does not carry them")
            pick += cand
    if len(pick) != sum(want.values()):
        raise RuntimeError(f"the seed pick is short: {pick}")
    dm = (c[:, closed] * 2.0) @ c[:, closed].conj().T
    dm += c[:, pick] @ c[:, pick].conj().T
    return dm.real.astype(np.float64)


def main() -> None:
    spec = campaign_spec("U14e")
    cfg = spec.atom_config()
    scratch = Path("elements/U14e/rungs/probe_tmp")
    scratch.mkdir(parents=True, exist_ok=True)
    for f in scratch.glob("*.chk"):        # a killed prior run's
        f.unlink()                         # scratch must not cache
    mol, mf = run_rohf(cfg, "elements/U14e/rungs/rohf.chk", ref=spec.ref)
    projs = lprojs(mol)
    print(f"head (trusted): E = {mf.e_tot:.8f}  open-block "
          f"{fmt(fingerprint(mol, mf, projs))}")

    for tag, closed, opens in MEMBERS:
        print(f"\n== member {tag}  closed={closed} open={opens} "
              f"spin={cfg.spin}")
        for acc in ACCELERATORS:
            for label, guess in ROHF_CENSUS_GUESSES:
                m = plain_arm(mol, closed, opens, acc, guess)
                report(f"{acc} {label}", mol, m, closed, opens, projs)
        # the seed arms: the configuration asked for explicitly (the
        # head's own MOs re-assigned; the hook stays, the seed owns
        # the start -- run_rohf's dm0= master), then the BARE seed
        # arm (no hook: plain aufbau maintains the seed) -- the
        # landscape measurement that says whether a stationary point
        # of the configuration exists at all, apart from the weight
        # hook's own preference
        want = {"l3": 4} if tag == "f4" else {"l3": 2, "l2": 2}
        dm = seed_density(mf, want)
        for acc in ACCELERATORS:
            _m, m = run_rohf(
                cfg, str(scratch / f"{tag}-{acc.replace('+', '-')}.chk"),
                ref=ReferenceSpec(occ_closed=closed, occ_open=opens,
                                  conv_accelerator=acc),
                dm0=dm)
            report(f"{acc} seed", mol, m, closed, opens, projs)
        _m, m = run_rohf(cfg, str(scratch / f"{tag}-bare.chk"),
                         ref=None, dm0=dm)
        print(f"  {'bare seed (aufbau)':>20}: E = {m.e_tot:.8f} "
              f"conv={bool(m.converged)}  | "
              f"{fmt(fingerprint(mol, m, projs))}", flush=True)
        if not m.converged:
            # the stall rescue (the twins' recipe's last leg): DIIS
            # held a neighborhood and stalled -- the Newton switch
            # converges the tail.  Two forms: plain Newton under
            # aufbau, and the hook re-attached over the stalled
            # density (its own occupied columns are the premium's
            # purest -- the ties break toward the stall's state)
            n = m.newton()
            n.conv_tol = 1e-12
            n.run(dm0=m.make_rdm1())
            print(f"  {'stall->newton (aufbau)':>20}: "
                  f"E = {n.e_tot:.8f} conv={bool(n.converged)}  | "
                  f"{fmt(fingerprint(mol, n, projs))}", flush=True)
            _m, n2 = run_rohf(
                cfg, str(scratch / f"{tag}-stall-hook.chk"),
                ref=ReferenceSpec(occ_closed=closed, occ_open=opens,
                                  conv_accelerator="newton"),
                dm0=m.make_rdm1())
            print(f"  {'stall->hook+newton':>20}: "
                  f"E = {n2.e_tot:.8f} conv={bool(n2.converged)} "
                  f"occ_ok={occupation_postcheck(mol, n2, closed, opens)}"
                  f"  | {fmt(fingerprint(mol, n2, projs))}", flush=True)

    for f in scratch.glob("*.chk"):
        f.unlink()
    scratch.rmdir()
    print("\ndone -- the declared basin per member: the deepest "
          "converged occ_ok arm carrying its configuration's "
          "fingerprint")


if __name__ == "__main__":
    main()
