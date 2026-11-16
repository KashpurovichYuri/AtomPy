"""The per-contraction dE-vs-anchor certificate table (ticket 37
Stage 3): for each distilled prior window (@r1/@s1/@s2 over the 3d
and 3d4s spans), the confined-valence mean-field cost against the
configuration's own trusted anchor -- Sc first, Ti mirrored.

Rows per twin:
  3d4s x {head, mixed} x {@r1, @s1, @s2}  (the full-valence hosts)
  3d   x {pure-d}      x {@r1, @s1, @s2}  (the d-manifold hosts)
  whole-space x {head, mixed, pure}       (the exactness legs)
+ the ion cross-check note (the twins' confined_core_orbitals
refusal, recorded).

Run: OMP_NUM_THREADS=16 PYTHONPATH=src python sandbox/table.py Sc-unc
The log lands in elements/<twin>/logs/subspace.log (committed).
"""

import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from atomterm.basis import series_window
from atomterm.cases import CAMPAIGN
from atomterm.system import build_mol
from subspace import ao_family, load_substrate, subspace_rohf

CORE_SHELLS = ("1s", "2s", "2p", "3s", "3p")


def main(key: str) -> None:
    cs = CAMPAIGN[key]
    cfg = cs.atom_config()
    mol = build_mol(cfg)
    s = mol.intor_symmetric("int1e_ovlp")
    fam = ao_family(mol, CORE_SHELLS, cs.shell_bands)
    head = load_substrate(f"elements/{key}/rungs/rohf.chk")
    members = {spec.tag: (spec, load_substrate(
        f"elements/{key}/rungs/series/{spec.tag}.chk"))
        for spec in cs.series[1:]}

    def spec_of(tag):
        if tag == "head":
            return (cs.ref.occ_closed, cs.ref.occ_open,
                    cs.ref.occ_e_ref, head)
        sp, mf = members[tag]
        return (sp.closed, sp.open, sp.e_ref, mf)

    # the window builders (the production paths; the seeds' mfs)
    mixed_tag = [t for t in members if t.endswith("s1")][0]
    pure_tag = [t for t in members if not t.endswith("s1")][0]
    all_mfs = [head] + [m for _, m in members.values()]

    def window(span_shells, tag, span_mfs):
        k = int(tag[1:])
        return series_window(mol, span_mfs, s, span_shells,
                             cs.prior_bands, {sh: k for sh in span_shells},
                             cs.ncore)[0]

    lines = [f"== {key}: the subspace-ROHF dE certificate "
             f"(2026-11-16, sandbox/table.py) ==",
             f"anchors: head {cs.ref.occ_e_ref:.9f}; " +
             "; ".join(f"{t} {members[t][0].e_ref:.9f}"
                       for t in (mixed_tag, pure_tag)),
             ""]
    print("\n".join(lines), flush=True)

    # the exactness legs (the whole space)
    n = mol.nao_nr()
    fam_rows = set(np.flatnonzero(np.any(np.abs(fam) > 1e-12, axis=1))
                   .tolist())
    rest = [i for i in range(n) if i not in fam_rows]
    eye = np.eye(n)
    lines.append("-- the whole-space exactness legs (the degenerate "
                 "limit: subspace == free ROHF) --")
    for tag in ("head", mixed_tag, pure_tag):
        closed, open_, anchor, mf = spec_of(tag)
        t0 = time.time()
        r = subspace_rohf(mf, fam, cs.ncore, eye[:, rest], closed, open_,
                          conv_accelerator="newton")
        d = (r["e_tot"] - anchor) * 1e3
        row = (f"whole      {tag:5s}  E {r['e_tot']:.9f}  "
               f"dE {d:+.6f} mEh  fill {r['report']['n_fill']}  "
               f"({time.time() - t0:.0f}s)")
        lines.append(row)
        print(row, flush=True)

    # the confined certificate rows
    lines.append("-- the confined legs (dE = E(subspace) - E(free "
                 "anchor), the window's mean-field cost) --")
    table = {}
    for span_shells, span_name, fills in (
            (("3d", "4s"), "3d4s", ("head", mixed_tag)),
            (("3d",), "3d", (pure_tag,))):
        for tag in fills:
            closed, open_, anchor, mf = spec_of(tag)
            win_mfs = [head] if False else all_mfs  # @r uses the head
            for wtag in ("r1", "s1", "s2"):
                if wtag == "r1":
                    win = window(span_shells, wtag, [head])
                else:
                    win = window(span_shells, wtag, all_mfs)
                t0 = time.time()
                r = subspace_rohf(mf, fam, cs.ncore, win, closed, open_,
                                  conv_accelerator="newton")
                d = (r["e_tot"] - anchor) * 1e3
                table[(span_name, tag, wtag)] = d
                row = (f"{span_name:5s} {wtag:3s} {tag:5s}  "
                       f"E {r['e_tot']:.9f}  dE {d:+.3f} mEh  "
                       f"capture {min(r['report']['capture']):.3f}  "
                       f"({time.time() - t0:.0f}s)")
                lines.append(row)
                print(row, flush=True)

    # the nesting note: within-@s family (the same eigenproblem, k
    # grows) vs cross-family (r1 vs s1 -- different eigenproblems,
    # nesting NOT expected; the spans' overlap is the identity axis)
    lines.append("-- the nesting read --")
    for span_name, tags in (("3d4s", ("head", mixed_tag)),
                            ("3d", (pure_tag,))):
        for tag in tags:
            d_s1 = table[(span_name, tag, "s1")]
            d_s2 = table[(span_name, tag, "s2")]
            d_r1 = table[(span_name, tag, "r1")]
            ok = "OK" if d_s2 <= d_s1 + 1e-6 else "VIOLATED"
            lines.append(f"{span_name:5s} {tag:5s}  dE(s2) - dE(s1) = "
                         f"{d_s2 - d_s1:+.3f} mEh  monotone: {ok};  "
                         f"dE(r1) = {d_r1:+.3f} (cross-family: the "
                         "head-block vs series-mean eigenproblems -- "
                         "nesting not expected)")
            print(lines[-1], flush=True)

    lines.append("-- the ion cross-check: confined_core_orbitals "
                 "REFUSES on the twins' family span (the closed-core "
                 "Roothaan floors at ~5e-10 past its own 100-cycle "
                 "budget -- the recorded declared class); the "
                 "cross-check runs on the carbon-class span "
                 "(sandbox/battery.py carbon: d = 0.000000 mEh) --")
    log = Path("elements") / key / "logs"
    log.mkdir(parents=True, exist_ok=True)
    (log / "subspace.log").write_text("\n".join(lines) + "\n")
    print(f"wrote {log / 'subspace.log'}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "Sc-unc")
