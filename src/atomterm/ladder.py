"""The orbital-treatment ladder: every rung through one analysis pass,
compared all-vs-all (ticket 08 on carbon, ticket 09 on Si).

Rungs, in order of increasing orbital freedom:

1. casci_raw    -- CASCI on the raw shell window, reference core folded as-is
                   (no spherical averaging): the deliberately symmetry-broken
                   reference. Fingerprint: [H, Lz] ~ eps and the multiplet
                   split {M_L=0} vs {+-1} (quadrupole k=2, m=0 core field;
                   the dipole k=1 is excluded by parity). In a sector that
                   is a complete single irrep L^2 = L(L+1)*Id identically,
                   so [H, L^2] carries no information there.
2. casci_avg    -- spherically averaged core: exact zeros. The scalar
                   projection of the polarized reference core.
3. casci_frac   -- average-of-configuration core (fractional-occupation SCF),
                   folded WITHOUT averaging: spherical at the source.
4. casscf_relax -- window span frozen; only core<->virtual rotations: adds
                   the core's polarization by the correlated active density.
5. casscf_radial -- window slots may rotate ONLY with same-(l, m) partners
                   (pure radial mixing, member-preserving basis + mask):
                   how much correlation is radially reachable WITHOUT any
                   angular freedom. NOTE exact closure additionally needs
                   the radial rotation applied identically to all m
                   partners; independent per-member radial mixing can leave
                   [H, Lz] at the member-radial-difference level (that is
                   the measurement).
6. casscf_virtmix -- core<->window exchange blocked; window may borrow
                   virtual character: the active pair's radial correlation,
                   isolated from identity-changing rotations.
7. casscf_shell -- free CASSCF from the raw-shell start.
8. casscf_pyscf -- free CASSCF from PySCF's default start (start dependence).
9. sa / sa_radial / sa_pyscf -- SA-CASSCF over the given weights: free from
                   the raw-shell start; constrained to radial rotations
                   (imposed symmetry -- does SA + the constraint drive BOTH
                   commutators to zero?); free from the canonical start.
10. *_pur rungs (optional, ``run_ladder(..., purified=True)``) -- the same
                   treatments on the PURIFIED core substrate
                   (:func:`purified_core_orbitals`: an exactly L-closed
                   core span): the fact-17 input-window baseline is gone,
                   so frozen-window rungs (relax, sa_relax) close at
                   machine precision on every atom, and free SA closes
                   exactly wherever the optimizer has a restoring force
                   (C/Si: yes; Sc: the flat member gauge leaves a ~1e-6
                   residual; the free-sa 4s<->3p swap attractor persists).

Adding a rung = appending one :class:`Rung` to :func:`_registry` (data, not
code): the runner executes it through :func:`stage_facts` and the report
iterates the registry, so a new rung can never silently vanish from the
report (that once happened to ``sa_pyscf``).

Divergence guard: every rung reports ``ortho_err`` = max|C^dag S C - I|.
Energies from rungs with ortho_err > 1e-6 are meaningless (the orbital
optimizer lost orthonormality -- observed on Si from raw-shell starts,
where reported energies plunge BELOW the variational bound). Such rungs
are flagged DIVERGED in the report and their numbers must not be read as
physics.

One-particle columns: ``|[P,Lz]|`` and ``|[P,Lp]|`` are the closure
defects of the window SPAN itself (operators.window_closure) -- no core
fold, no CI. They attribute each rung's many-body commutators to the
one-particle space: ``|[P,Lz]|`` is the part [H, Lz] can see (blind to
m=0-member and equal-+-pair tilts), ``|[P,Lp]|`` the part only L^2 sees.
Both zero = the complete closure certificate of the span (fact 19).

max|dE| dictionary (our reconstruction vs PySCF's own energies):
  0.0          -- identical folds (raw/frac rungs), or a spherical core
                   (SA: the ensemble density is a rotational scalar);
  ~1e-5..5e-5  -- the energy footprint of a POLARIZED core that we average
                   and PySCF uses raw (avg, optimized rungs). A polarization
                   detector, not an error.

Sector subtlety: in a maximal-Ms sector that is a complete single irrep
(carbon/Si 3-det Ms=1), [H, L^2] = 0 for ANY H -- only [H, Lz] and the
splitting are fingerprints there. In sectors mixing several terms both
commutators are informative.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np

from atomterm.system import AtomConfig
from atomterm.drivers import (
    build_full_mo,
    constrained_casscf,
    core_virtual_only_mask,
    purified_core_orbitals,
    radial_basis,
    radial_mask,
    run_casci,
    run_casscf,
    run_casscf_core_relax,
    run_casscf_pyscf_start,
    run_casscf_virtual_mix,
    run_frac_occ_scf,
    run_rohf,
    run_sa_casscf,
    select_valence_ao,
)
from atomterm.stages import StageFacts, stage_facts

__all__ = ["Rung", "run_carbon_ladder", "run_ladder"]


@dataclass(frozen=True)
class Rung:
    """One ladder rung: a stage to run plus how to fold-and-report it.

    tag : rung name (report row, chkfile suffix, payload key).
    kind : reporting kind ("casci" / "casscf") -- CASCI rungs also report
        ``core_anisotropy`` of their reference core.
    run : executes the stage, returns a CASResult (chkfile-cached).
    window : active columns to fold on -- raw shell AOs for CASCI rungs,
        the converged ``mo_coeff`` window for optimized rungs.
    fold_mf : reference mf object used in the fold (the fractional rung
        folds on its own mf_f; integrals are geometry-defined anyway).
    spherical_core : fold with the averaged core (True) or raw (False).
    anisotropy_ref : ``(mf_ref, c_core_or_None)`` for CASCI reference rungs;
        ``None`` skips the anisotropy column.
    """

    tag: str
    kind: str
    run: Callable[[], object]
    window: Callable[[object], np.ndarray]
    fold_mf: Callable[[], object]
    spherical_core: bool = True
    anisotropy_ref: Callable[[], tuple] | None = None


def _rung_payload(rung: Rung, cas, mol, ncore: int, na: int, nb: int) -> dict:
    """Execute one rung's analysis pass: run -> fold -> facts -> payload."""
    facts: StageFacts = stage_facts(
        mol, rung.fold_mf(), cas, rung.window(cas), ncore, na, nb,
        tag=rung.tag, spherical_core=rung.spherical_core,
        anisotropy_ref=(rung.anisotropy_ref() if rung.anisotropy_ref else None),
    )
    return facts.payload(rung.kind)


def _registry(mf, mf_f, ncore: int, ncas: int, nelecas: tuple[int, int],
              act_idx: list[int], ndet: int, sa_weights, chk,
              c_pur: np.ndarray = None) -> list[Rung]:
    """The ladder as data. ``chk`` maps a rung tag to its chkfile.

    ``c_pur`` (optional): purified core columns
    (:func:`purified_core_orbitals`) -- appends the ``*_pur`` rungs, the
    same orbital treatments on the exactly-closed substrate (the fact-17
    input-window baseline removed; ticket 13)."""
    n_open = len(act_idx)
    c_shell = lambda cas: np.eye(mf.mo_coeff.shape[0])[:, act_idx]  # noqa: E731
    c_win = lambda cas: cas.mo_coeff[:, ncore : ncore + ncas]  # noqa: E731
    mo_of = lambda mf_ref, c_core=None: build_full_mo(  # noqa: E731
        mf_ref, list(range(ncore)), act_idx, c_core=c_core)
    nroots = min(9, ndet)

    rungs = [
        Rung(
            "casci_raw", "casci",
            lambda: run_casci(mf, ncas, nelecas, mo_of(mf), chk("casci_raw"),
                              nroots=nroots, name="casci_raw"),
            c_shell, lambda: mf, spherical_core=False,
            anisotropy_ref=lambda: (mf, None),
        ),
        Rung(
            "casci_avg", "casci",
            lambda: run_casci(mf, ncas, nelecas, mo_of(mf), chk("casci_avg"),
                              nroots=nroots, name="casci_avg"),
            c_shell, lambda: mf,
            anisotropy_ref=lambda: (mf, None),
        ),
    ]
    if mf_f is not None:
        c_frac = lambda: mf_f.mo_coeff[:, :ncore]  # noqa: E731
        rungs.append(
            Rung(
                "casci_frac", "casci",
                lambda: run_casci(mf_f, ncas, nelecas, mo_of(mf_f, c_frac()),
                                  chk("casci_frac"), nroots=nroots,
                                  name="casci_frac"),
                c_shell, lambda: mf_f,
                anisotropy_ref=lambda: (mf_f, c_frac()),
            )
        )

    mo_shell = mo_of(mf)
    # member-labeled basis + same-(l, member) rotation mask: the radial rungs
    mo_rad, lab_rad = radial_basis(mf, list(range(ncore)), act_idx)
    mask_rad = radial_mask(lab_rad, ncore, ncas)
    rungs += [
        Rung(
            "casscf_relax", "casscf",
            lambda: run_casscf_core_relax(
                mf, ncas, nelecas, mo_shell, chk("casscf_relax")),
            c_win, lambda: mf,
        ),
        Rung(
            "casscf_radial", "casscf",
            lambda: constrained_casscf(
                mf, ncas, nelecas, chk("casscf_radial"), mo_full=mo_rad,
                mask=mask_rad, name="casscf_radial"),
            c_win, lambda: mf,
        ),
        Rung(
            "casscf_virtmix", "casscf",
            lambda: run_casscf_virtual_mix(
                mf, ncas, nelecas, mo_shell, chk("casscf_virtmix")),
            c_win, lambda: mf,
        ),
        Rung(
            "casscf_shell", "casscf",
            lambda: run_casscf(
                mf, ncas, nelecas, mo_shell, chk("casscf_shell")),
            c_win, lambda: mf,
        ),
        Rung(
            "casscf_pyscf", "casscf",
            lambda: run_casscf_pyscf_start(
                mf, ncas, nelecas, chk("casscf_pyscf")),
            c_win, lambda: mf,
        ),
    ]
    if sa_weights:
        rungs += [
            Rung(
                "sa", "casscf",
                lambda: run_sa_casscf(
                    mf, ncas, nelecas, mo_shell, list(sa_weights), chk("sa")),
                c_win, lambda: mf,
            ),
            # imposed-symmetry SA: same weights, rotations restricted to
            # same-(l, member) radial mixing -- the constrained twin of "sa"
            Rung(
                "sa_radial", "casscf",
                lambda: run_sa_casscf(
                    mf, ncas, nelecas, mo_rad, list(sa_weights),
                    chk("sa_radial"), mask=mask_rad),
                c_win, lambda: mf,
            ),
            # canonical-start variant: isolates the optimizer's start
            # dependence (Si's raw-shell start diverges; the canonical
            # one converges)
            Rung(
                "sa_pyscf", "casscf",
                lambda: run_sa_casscf(
                    mf, ncas, nelecas, mf.mo_coeff, list(sa_weights),
                    chk("sa_pyscf")),
                c_win, lambda: mf,
            ),
        ]
    if c_pur is not None:
        # the purified-substrate twins (ticket 13): same treatments, core
        # span exactly L-closed -> input windows exactly closed. The CASCI
        # rung folds the window PySCF actually diagonalized (c_win, not
        # the raw AOs) AND its purified core density is m-symmetric by
        # construction, so raw fold == averaged fold: max|dE| drops to
        # machine precision (the fact-14 window-definition offset and the
        # anisotropy detector both vanish on this substrate).
        mo_shell_p = mo_of(mf, c_pur)
        mo_rad_p, lab_rad_p = radial_basis(
            mf, list(range(ncore)), act_idx, c_core=c_pur)
        mask_rad_p = radial_mask(lab_rad_p, ncore, ncas)
        rungs += [
            Rung(
                "casci_avg_pur", "casci",
                lambda: run_casci(mf, ncas, nelecas, mo_shell_p,
                                  chk("casci_avg_pur"), nroots=nroots,
                                  name="casci_avg_pur"),
                c_win, lambda: mf,
                anisotropy_ref=lambda: (mf, c_pur),
            ),
            Rung(
                "casscf_relax_pur", "casscf",
                lambda: run_casscf_core_relax(
                    mf, ncas, nelecas, mo_shell_p, chk("casscf_relax_pur")),
                c_win, lambda: mf,
            ),
            Rung(
                "casscf_radial_pur", "casscf",
                lambda: constrained_casscf(
                    mf, ncas, nelecas, chk("casscf_radial_pur"),
                    mo_full=mo_rad_p, mask=mask_rad_p,
                    name="casscf_radial_pur"),
                c_win, lambda: mf,
            ),
        ]
        if sa_weights:
            rungs += [
                # THE exact-symmetry SA rung: window frozen (closure exact
                # by construction), core relaxed under the ensemble
                Rung(
                    "sa_relax_pur", "casscf",
                    lambda: run_sa_casscf(
                        mf, ncas, nelecas, mo_shell_p, list(sa_weights),
                        chk("sa_relax_pur"), mask=core_virtual_only_mask),
                    c_win, lambda: mf,
                ),
                Rung(
                    "sa_pur", "casscf",
                    lambda: run_sa_casscf(
                        mf, ncas, nelecas, mo_shell_p, list(sa_weights),
                        chk("sa_pur")),
                    c_win, lambda: mf,
                ),
                Rung(
                    "sa_radial_pur", "casscf",
                    lambda: run_sa_casscf(
                        mf, ncas, nelecas, mo_rad_p, list(sa_weights),
                        chk("sa_radial_pur"), mask=mask_rad_p),
                    c_win, lambda: mf,
                ),
            ]
    return rungs


def run_ladder(cfg: AtomConfig, shell_patterns: list[str], nelec_act: int,
               workdir: str = "results", sa_weights: list[float] | None = None,
               frac_nelec_open: int | None = None, frac_n_open: int | None = None,
               rung_tags: list[str] | None = None,
               purified: bool = False) -> dict:
    """Run ladder rungs for one atom's primary Ms sector, print and persist
    the all-vs-all comparison. Returns ``{rung_tag: payload}``.

    rung_tags : optional subset of rung tags to run (e.g. window studies
        that only need the CASCI + free-CASSCF rungs); ``None`` = all.
    purified : append the ``*_pur`` rungs (the same treatments on the
        exactly-closed purified core substrate, ticket 13 / fact 17).
    """
    out = Path(workdir)
    out.mkdir(parents=True, exist_ok=True)
    label = cfg.label

    mol, mf = run_rohf(cfg, str(out / f"{label}_rohf.chk"))
    act_idx = select_valence_ao(mol, shell_patterns)
    mol_f = mf_f = None
    if frac_nelec_open is not None:
        # open_idx: select the fractional block by CHARACTER (the active
        # shell AOs). Energy-ordered filling is fine for C/Si but runs away
        # on TMs where 3d/4p quasi-degenerate (Sc: p-containing ensemble).
        mol_f, mf_f = run_frac_occ_scf(
            cfg, str(out / f"{label}_fracocc.chk"), frac_nelec_open, frac_n_open,
            open_idx=act_idx)
    ncas = len(act_idx)
    spin = cfg.spin
    na, nb = (nelec_act + spin) // 2, (nelec_act - spin) // 2
    ncore = (mol.nelectron - nelec_act) // 2
    chk = lambda name: str(out / f"{label}_{name}.chk")  # noqa: E731
    ndet = math.comb(ncas, na) * math.comb(ncas, nb)

    c_pur = None
    if purified:
        c_pur, _ = purified_core_orbitals(mol, mf, ncore)

    registry = _registry(mf, mf_f, ncore, ncas, (na, nb), act_idx, ndet,
                         sa_weights, chk, c_pur=c_pur)
    if rung_tags is not None:
        wanted = set(rung_tags)
        unknown = wanted - {r.tag for r in registry}
        if unknown:
            raise ValueError(
                f"unknown rung tags {sorted(unknown)}; available: "
                f"{[r.tag for r in registry]}")
        registry = [r for r in registry if r.tag in wanted]

    rungs: dict[str, dict] = {}
    for rung in registry:
        rungs[rung.tag] = _rung_payload(rung, rung.run(), mol, ncore, na, nb)

    _report(rungs, out, label, mf)
    return rungs


def run_carbon_ladder(cfg: AtomConfig | None = None, workdir: str = "results") -> dict:
    """The carbon ladder (ticket 08): 2p^2 through every rung."""
    return run_ladder(cfg or AtomConfig(), ["2p"], 2, workdir=workdir,
                      sa_weights=[1 / 3] * 3, frac_nelec_open=2, frac_n_open=3)


def _report(rungs: dict, out: Path, label: str, mf) -> None:
    """Print and persist the per-rung comparison table + all-vs-all matrix."""
    # the report follows the registry (= insertion) order, so it can never
    # drop a rung that ran
    tags = list(rungs)
    ref = rungs.get("casci_avg")

    lines = [f"# Orbital-treatment ladder: {label}", "",
             f"ROHF reference energy: {mf.e_tot:.8f} Eh", "",
             f"{'rung':<14} {'E_ground':>16} {'dE_vs_avg':>10} {'<L^2>':>9} "
             f"{'|[H,Lz]|':>10} {'|[H,L2]|':>10} {'|[P,Lz]|':>9} {'|[P,Lp]|':>9} "
             f"{'win_l':>6} {'aniso':>9} {'max|dE|':>9} {'ortho':>9}  ground config"]
    for tag in tags:
        r = rungs[tag]
        d_e = r["e_ground"] - ref["e_ground"] if ref else 0.0
        purity = max(r["window"].values())
        aniso = r.get("core_anisotropy")
        flag = " DIVERGED" if r["diverged"] else ""
        lines.append(
            f"{tag:<14} {r['e_ground']:>16.8f} {d_e:>10.3e} {r['ground_l2']:>9.5f} "
            f"{r['comm_norm_h_lz']:>10.2e} {r['comm_norm_h_l2']:>10.2e} "
            f"{r['win_lz']:>9.1e} {r['win_lperp']:>9.1e} "
            f"{purity:>6.3f} {aniso if aniso is not None else float('nan'):>9.1e} "
            f"{r['max_de_vs_pyscf']:>9.1e} {r['ortho_err']:>9.1e}  "
            f"{r['configurations'][0]}{flag}"
        )

    lines += ["", "pairwise ground-energy differences (E_i - E_j), Eh:",
              "            " + " ".join(f"{t[:11]:>12}" for t in tags)]
    for ti in tags:
        row = [f"{rungs[ti]['e_ground'] - rungs[tj]['e_ground']:>12.3e}"
               for tj in tags]
        lines.append(f"{ti:<12}" + " ".join(row))

    lines += ["", "configurations per root:"]
    for tag in tags:
        for i, cfg_str in enumerate(rungs[tag]["configurations"]):
            lines.append(f"  {tag:<14} root {i}: {cfg_str}")

    text = "\n".join(lines)
    print(text)
    (out / f"{label}_ladder.md").write_text(text + "\n", encoding="utf-8")
    (out / f"{label}_ladder.json").write_text(
        json.dumps(rungs, indent=2), encoding="utf-8")
