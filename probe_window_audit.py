"""Probe: the window-composition audit of the optimizer checkpoints
(scratch; the log is the record).  The l-space weight (shell_weights,
the certified partition) aggregated over the declared bands' angular
classes: the span's on-band weight and the escape in column units.
The per-band |c|^2 profile = the qualitative radial record.
Usage: python probe_window_audit.py <addr>:<ncore>:<bands-comma> ...
"""
import re
import sys
from pathlib import Path

import numpy as np
from pyscf import lib

sys.path.insert(0, str(Path(__file__).resolve().parents[0] / "src"))

from atomterm.cases import campaign_spec  # noqa: E402
from atomterm.windows import select_shell_ao, shell_weights  # noqa: E402


def band_ao_groups(mol):
    groups = {}
    for mu, entry in enumerate(mol.ao_labels(fmt=False)):
        m = re.match(r"(\d+[spdfgh])", entry[2])
        if m is None:
            raise ValueError(f"AO label {entry[2]!r} has no band token")
        groups.setdefault(m.group(1), []).append(mu)
    return groups


def main(argv):
    spec = campaign_spec("U14e")
    from atomterm.reference import run_rohf
    mol, mf = run_rohf(spec.atom_config(),
                       "elements/U14e/rungs/rohf.chk", ref=spec.ref)
    groups = band_ao_groups(mol)
    s = mf.get_ovlp()

    for arg in argv:
        addr, ncore, bandstr = arg.split(":")
        ncore = int(ncore)
        declared = bandstr.split(",")
        chk = f"elements/U14e/rungs/{addr}.chk"
        mo = lib.chkfile.load(chk, f"atomterm/{addr}/mo_coeff")
        ncas = int(lib.chkfile.load(chk, f"atomterm/{addr}/ncas"))
        act = np.asarray(mo)[:, ncore:ncore + ncas]
        names, w = shell_weights(mol, act, s)
        lcol = {nm: k for k, nm in enumerate(names)}
        onband_l = sorted({b[-1] for b in declared})
        i_on = [lcol[l] for l in onband_l]
        w_on = w[:, i_on].sum(axis=1)
        print(f"== {addr}  (ncore={ncore}, ncas={ncas}, declared "
              f"{bandstr}; on-band l = {''.join(onband_l)}) ==", flush=True)
        prof = {b: np.sum(np.abs(act[rows, :]) ** 2, axis=0).real
                for b, rows in groups.items()}
        for j in range(ncas):
            top = sorted(groups, key=lambda b: -prof[b][j])[:2]
            p = "  ".join(f"{b}:{prof[b][j]:.3f}" for b in top)
            ls = "  ".join(f"{nm}:{w[j, lcol[nm]]:.3f}" for nm in names
                           if w[j, lcol[nm]] > 0.01)
            print(f"  col[{j:2d}]  on-l {w_on[j]:.3f}   [{ls}]   {p}",
                  flush=True)
        frac = float(np.sum(w_on))
        print(f"  SPAN: on-band l-weight {frac:.3f} of {ncas} -> escape "
              f"{ncas - frac:.2f} columns outside the declared "
              f"angular classes", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
