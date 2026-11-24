"""The channel-map probe: which undeclared primitive bands carry the
reference's density and the low virtuals -- the pre-CAS viability map
for the correlation-channel exploration (prints only).

THREE SIGNALS per primitive band (generic instruments, any system):
  1. the occupied POPULATION (the ROHF occupied MOs' band weights) --
     the static-density view;
  2. the excitation cost -- the min virtual MO energy among the MOs
     dominated by the band -- the intruder/cheapest-channel view;
  3. the device's radial spectrum over each shell FAMILY (s/p/d/f,
     the occupied-signal Gram, the m-summed device quantity) with the
     per-band weight of each radial -- the channel map proper: which
     primitives the significant radials live on.

Everything read from the trusted rung (the C1 ROHF).  No CAS runs.

Run from the repo root:
    OMP_NUM_THREADS=8 python scripts/probes/probe_u14e_channel_map.py \
        | tee logs/probe_u14e_channel_map.log
"""

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "u14e"))

import numpy as np  # noqa: E402

from _vv_common import (band_chains_by_suffix, band_population,  # noqa: E402
                        family_spectrum, load_reference)

DECLARED = {"5f": ("9f", "10f", "11f"), "6d": ("12d", "13d", "14d"),
            "6p": ("13p", "14p", "15p"), "7s": ("18s", "19s", "20s")}
NCORE = 5   # the folded-core pairs: the valence density cut (the
            # device's ncol convention at the (4,0) partition)


def main() -> int:
    t00 = time.time()
    mol, mf = load_reference()
    occ = np.asarray(mf.mo_occ) > 0
    e_all = np.asarray(mf.mo_energy)

    # ---- the bands and the MO band decomposition ------------------------
    bands, types = band_chains_by_suffix(mol)
    order, frac = band_population(mol, mf, bands, weights=mf.mo_occ)
    dom = frac.argmax(axis=1)          # each MO's dominant band index
    print(f"  [bands] {len(order)} primitive bands", flush=True)

    # ---- per-primitive: the occupied population and the min virtual ε
    occ_pop = np.zeros(len(order))
    min_virt = {}
    virt_dom = {}
    for j in range(frac.shape[0]):
        k = dom[j]
        occ_pop[k] += (mf.mo_occ[j] if occ[j] else 0.0)
        if not occ[j]:
            virt_dom.setdefault(k, []).append((e_all[j], j))
    for k, lst in virt_dom.items():
        min_virt[k] = min(e for e, _ in lst)

    print("\n===== the per-primitive map (sorted by the occupied "
          "population) =====", flush=True)
    print("  band       occ_e   min_virtual_Eh   declared?", flush=True)
    decl_of = {b: shell for shell, bs in DECLARED.items() for b in bs}
    for k in np.argsort(-occ_pop):
        b = order[k]
        mv = min_virt.get(k)
        decl = decl_of.get(b, "undeclared")
        print(f"  {b:6s}  {occ_pop[k]:7.3f}   "
              f"{('%.4f' % mv) if mv is not None else '—':>10s}   "
              f"{decl}", flush=True)

    # ---- the Lowdin MO table: the occupied + the first 25 virtuals ----
    print("\n===== the MO table (the occupied + the first 25 "
          "virtuals; the top-2 bands per MO) =====", flush=True)
    virt_sorted = sorted(np.flatnonzero(~occ), key=lambda j: e_all[j])
    show = list(np.flatnonzero(occ)) + list(virt_sorted[:25])
    for j in show:
        fr = frac[j]
        top = np.argsort(-fr)[:2]
        e = e_all[j]
        print(f"  MO {j:3d} {'occ' if occ[j] else 'vir'} eps "
              f"{e:+.4f}: "
              + "  ".join(f"{order[t]} {fr[t]:.2f}" for t in top),
              flush=True)

    # ---- stage 3: the family radial spectra -----------------------------
    print("\n===== the family radial spectra (the occupied-signal "
          "Gram, m-summed; the per-band weights of each radial) =====",
          flush=True)
    for shell, fam_bands in (("s", [b for b in order if b.endswith("s")]),
                             ("p", [b for b in order if b.endswith("p")]),
                             ("d", [b for b in order if b.endswith("d")]),
                             ("f", [b for b in order if b.endswith("f")])):
        lam, (ubands, wb) = family_spectrum(mol, mf, shell, fam_bands,
                                            NCORE, types)
        print(f"  [{shell}] the family slots' spectrum (leading 6): "
              + ", ".join(f"{x:.4f}" for x in lam[:6])
              + f"  |  the bands: {ubands}", flush=True)
        for e in range(min(3, len(lam))):
            w = wb[e]
            top = np.argsort(-w)[:4]
            print(f"    radial {e + 1} (Eig {lam[e]:.4f}): "
                  + "  ".join(f"{ubands[t]} {w[t]:.2f}" for t in top
                              if w[t] > 0.01), flush=True)

    # ---- the intruder screen --------------------------------------------
    print("\n===== the intruder screen (the channels hosting the "
          "lowest virtuals) =====", flush=True)
    e_min_all = min(min_virt.values()) if min_virt else float("nan")
    for k in np.argsort(np.array([min_virt.get(k, np.inf)
                                  for k in range(len(order))]))[:5]:
        b = order[k]
        mv = min_virt.get(k)
        if mv is None:
            continue
        print(f"  {b:6s}: the min virtual {mv:+.4f} Eh "
              f"(the occupied pop {occ_pop[k]:.3f}; declared: "
              f"{decl_of.get(b, 'undeclared')})", flush=True)

    print(f"\n({time.time() - t00:.0f} s total)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
