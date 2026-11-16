"""The subspace-ROHF bench (ticket 37 Stage 3).

The ``run_folded`` mirror INVERTED: the VALENCE is confined to a
distilled prior window (``@r1``/``@s1``/``@s2`` columns) while the CORE
stays exact -- its full declared AO family, free to relax (the cusp
lesson: the core's energy is band-complete, only the valence pays the
distillation).  The certificate: the per-contraction mean-field cost
dE = E(subspace) - E(free anchor), measured per window and per
configuration.

Mechanics: PySCF's ROHF kernel runs UNMODIFIED on the S-orthonormal
subspace ``[core family Loewdin'd | window projected against it]``
(``get_ovlp = I``, projected ``h``, the ``veff`` shim translating the
subspace density to the AO field and projecting the result back) --
the same shims :func:`atomterm.reference.run_folded` uses, so the
shared tiny helpers are imported from there.  Everything (core +
valence) sits INSIDE the SCF, so the kernel's own ``e_tot`` is the
honest restricted variational value: no custom bookkeeping.  The
filling is block-aware MOM: the core pairs among the family block, the
valence closed pairs and opens among the window block, every pick
anchored on the projected trusted reference (the seed declares the
configuration; the window block is too small for aufbau to be safe --
a d column dipping below 4s would silently flip it).  The seedless ion
leg (battery 7's closed-shell cross-check) fills the family block by
aufbau instead -- :func:`atomterm.core.confined_core_orbitals`'s own
rule, the two constructions meeting through different kernels.

Sandbox code, imported by nothing; integration lands this as
``reference.subspace_rohf`` (git history keeps the bench).
"""

import numpy as np
from pyscf import gto, scf

from atomterm.reference import _closed_field, _spin_summed
from atomterm.windows import lowdin, pattern_projector, realized_patterns, \
    select_valence_ao


def load_substrate(chk: str):
    """The trusted reference chk as a substrate mf (cache-trusted:
    converged columns, the real mol, the kernel accessors)."""
    mol, rec = scf.chkfile.load_scf(chk)
    mf = scf.ROHF(mol)
    mf.__dict__.update(rec)
    mf.converged = True
    return mf


def ao_family(mol, shells: tuple[str, ...], layer: dict) -> np.ndarray:
    """The shells' full AO family as raw AO columns -- the declared
    band realization (``layer``), the cusp-honest core space before
    orthonormalization."""
    rows = select_valence_ao(mol, realized_patterns(shells, layer))
    return np.eye(mol.nao_nr())[:, rows]


def build_subspace(mf, c_fam_raw: np.ndarray,
                   c_win: np.ndarray) -> tuple[np.ndarray, int]:
    """The S-orthonormal subspace ``[family Loewdin'd | window
    projected against it]``: the family orthonormalized, the window
    columns family-orthogonalized and re-orthonormalized.  Returns
    ``(Q, n_fam)``."""
    s = mf.get_ovlp()
    c_fam = lowdin(c_fam_raw, s)
    if c_win.shape[1] == 0:
        return c_fam, c_fam.shape[1]
    w = c_win - c_fam @ (c_fam.conj().T @ s @ c_win)
    w = lowdin(w, s)
    return np.hstack([c_fam, w]), c_fam.shape[1]


def subspace_rohf(mf, c_fam_raw: np.ndarray, n_pairs: int, c_win: np.ndarray,
                  occ_closed: tuple[str, ...], occ_open: tuple[str, ...],
                  conv_accelerator: str = "diis"):
    """One confined mean-field solve: the atom's full ROHF restricted
    to ``[core family | window]``.

    Parameters
    ----------
    mf : the trusted substrate (the seed's occupied columns carry the
        configuration; ``mol``, ``get_hcore``, ``get_ovlp``,
        ``get_veff`` for the kernel).
    c_fam_raw : (nao, M) the core family's raw AO columns (declared
        bands; Loewdin'd inside).
    n_pairs : the core's chemistry (pairs).
    c_win : (nao, W) the window columns (S-orthonormal from
        :func:`atomterm.basis.series_window`, or raw AO columns).
    occ_closed, occ_open : the filling gates' character patterns (the
        honest target-side l-vocabulary).

    Returns
    -------
    dict
        ``e_tot`` (the kernel's own restricted variational total),
        ``converged``, ``c_ao`` (the converged MO columns back in AO
        coordinates), ``mo_occ``, and the disclosure: family/window
        widths, per-seed-column capture, the block splits.

    Raises
    ------
    RuntimeError
        Non-convergence, a lost filling count, or a seed column that
        does not land in its declared block (capture below 0.5 -- the
        subspace cannot host the trusted configuration).
    """
    mol = mf.mol
    n_open = mol.spin
    n_closed_val = (mol.nelectron - 2 * n_pairs - mol.spin) // 2
    if n_closed_val < 0:
        raise ValueError(
            f"{n_pairs} core pairs leave {mol.nelectron - 2 * n_pairs} "
            f"valence electrons but the sector wants {mol.spin} open")
    Q, n_fam = build_subspace(mf, c_fam_raw, c_win)
    n_win = Q.shape[1] - n_fam
    if n_fam < n_pairs:
        raise ValueError(
            f"the core family holds {n_fam} columns; the core's "
            f"{n_pairs} pairs need more -- the family declaration is "
            "too narrow")
    if n_win < 2 * n_closed_val + n_open:
        raise ValueError(
            f"the window block holds {n_win} columns; the valence "
            f"filling wants {2 * n_closed_val + n_open}")

    s = mf.get_ovlp()
    # the seed: the trusted occupied columns projected into the
    # subspace (block-disjoint families land cleanly; the per-column
    # capture is the disclosure).  The closed/open anchor split
    # follows the trusted OCCUPATIONS (canonical column order ties
    # the 4s pair against the d open by energy -- the values carry
    # the configuration)
    has_seed = (getattr(mf, "mo_coeff", None) is not None
                and getattr(mf, "mo_occ", None) is not None)
    if has_seed:
        occ_seed = np.asarray(mf.mo_occ)
        occ_rows = np.flatnonzero(occ_seed > 0)
        occ_seed = occ_seed[occ_rows]
        c_seed = Q.conj().T @ s @ mf.mo_coeff[:, occ_rows]
        capture = np.linalg.norm(c_seed, axis=0) ** 2
        c_seed = c_seed / np.linalg.norm(c_seed, axis=0)
        if capture.min() < 0.5:
            raise RuntimeError(
                f"a seed column captures only {capture.min():.3f} of "
                "itself in the subspace -- the subspace cannot host "
                "the trusted configuration")
        # the seed density in SPIN-SEPARATED form: a summed dm would be
        # split (d/2, d/2) by the kernel's first-cycle handling -- the
        # open shell's spin asymmetry (the whole point of ROHF) lost,
        # the first Fock garbage, the filling off its basin from cycle
        # one (measured: the carbon whole-space leg converged into an
        # excited configuration, the Sc leg ran away)
        occ_a = np.minimum(occ_seed, 1.0)
        occ_b = occ_seed - occ_a
        dm0 = np.stack([(c_seed * occ_a) @ c_seed.conj().T,
                        (c_seed * occ_b) @ c_seed.conj().T])
    else:
        capture = None
        dm0 = None

    # core-block seed split: the columns whose projection lives in the
    # family block carry the core's pairs; for the seedless ion leg
    # everything fills by aufbau
    p_closed = pattern_projector(mol, occ_closed)
    p_open = pattern_projector(mol, occ_open)

    hc_sub = Q.conj().T @ mf.get_hcore() @ Q
    eye = np.eye(Q.shape[1])
    if not has_seed:
        # the ion leg's aufbau start on the subspace h, spin-separated
        # (a closed shell: both spins carry the same density)
        _ev, U = np.linalg.eigh(hc_sub)
        occ0 = np.zeros(Q.shape[1])
        occ0[:n_pairs] = 2.0
        dm = (U * occ0) @ U.conj().T
        dm0 = np.stack([dm * 0.5, dm * 0.5])

    mol_scf = mol
    if getattr(mol, "symmetry", False):
        mol_scf = gto.M(atom=mol.atom, basis=mol.basis, spin=mol.spin,
                        ecp=mol.ecp, charge=mol.charge, verbose=0)

    fm = scf.ROHF(mol_scf)
    fm.conv_tol = 1e-12
    fm.verbose = 0
    fm.get_hcore = lambda *args, **kw: hc_sub
    fm.get_ovlp = lambda *args, **kw: eye

    prev_cl = [None]
    prev_op = [None]
    if has_seed:
        # the MOM anchors: the seed's closed columns (core pairs AND
        # the valence closed -- both fill 2.0, the distinction carries
        # no filling weight) and its open columns
        prev_cl[0] = c_seed[:, occ_seed == 2.0]
        prev_op[0] = c_seed[:, occ_seed == 1.0]

    # the initial density (the seed pair or the ion's aufbau start) --
    # set for BOTH branches: the kernel's fallback guess builds in the
    # REAL mol's AO space, wrong coordinates for the shimmed kernel
    def _init_guess(mol_=None, key=None, **kw):
        return dm0

    fm.get_init_guess = _init_guess

    def _w_mask(mo_coeff, P):
        """Each MO column's span weight: the projector through the
        subspace (the FULL coefficient columns -- slicing rows would
        project fragments, not orbitals)."""
        c_ao = Q @ mo_coeff
        return np.einsum("ai,ab,bi->i", c_ao.conj(), P, c_ao).real

    n_cl_tot = n_pairs + n_closed_val

    def _mom_assign(prev_cols, mo_coeff, allowed, n_pick, mo_energy):
        """The certified MOM assignment (_mom_occ's form): each seed
        slot, in seed order, takes its best-overlap REMAINING column
        -- a per-slot contest the physical columns win by
        self-overlap (a summed-anchor score lets the metric-amplified
        garbage, whose directions live inside the closed spans, win
        the aggregate -- measured: the head whole-space leg diverged
        to +1e7 under a summed form).  Shortfalls fill by aufbau."""
        picked = []
        used = set()
        for j in range(prev_cols.shape[1]):
            if len(picked) >= n_pick:
                break
            free = [i for i in allowed if i not in used]
            if not free:
                break
            ovl = np.abs(prev_cols[:, j].conj() @ mo_coeff)
            picked.append(max(free, key=lambda i: ovl[i]))
            used.add(picked[-1])
        while len(picked) < n_pick:
            free = [i for i in allowed if i not in used]
            if not free:
                break
            picked.append(min(free, key=lambda i: mo_energy[i]))
            used.add(picked[-1])
        return picked

    # the last projected pair (the clean-Fock source for the fills)
    vhf_cell = [None]

    def _veff(mol_=None, dm=None, dm_last=0, vhf_last=None, hermi=1):
        # the density translates to the AO field and the TRUE
        # per-spin fields project back: the pair's exchange follows
        # EACH spin density (a (v, v) pair built from the summed
        # density carries the closed-shell exchange formula -- correct
        # only for closed shells, 0.3-Eh-class wrong on the opens)
        d = np.asarray(dm)
        if d.ndim == 3:
            d_pair = np.stack([Q @ di @ Q.conj().T for di in d])
            v = mf.get_veff(mol, d_pair)
            vp = np.stack([Q.conj().T @ v[0] @ Q,
                           Q.conj().T @ v[1] @ Q])
        else:
            d_ao = Q @ d @ Q.conj().T
            vv = Q.conj().T @ _closed_field(mf, d_ao) @ Q
            vp = np.array((vv, vv))
        vhf_cell[0] = vp
        return vp

    fm.get_veff = _veff

    def _jk(mol_=None, dm=None, hermi=1, with_j=True, with_k=True):
        """The J/K translator the Newton response functions contract
        through: subspace densities in, subspace J/K out (the AO
        build happens on the translated density, the fields project
        back)."""
        d = np.asarray(dm)
        if d.ndim == 3:
            jk = [mf.get_jk(mol, Q @ di @ Q.conj().T, hermi=hermi,
                            with_j=with_j, with_k=with_k) for di in d]
            vj = np.stack([Q.conj().T @ j @ Q for j, _ in jk])
            vk = np.stack([Q.conj().T @ k @ Q for _, k in jk])
            return vj, vk
        d_ao = Q @ d @ Q.conj().T
        vj, vk = mf.get_jk(mol, d_ao, hermi=hermi, with_j=with_j,
                           with_k=with_k)
        return Q.conj().T @ vj @ Q, Q.conj().T @ vk @ Q

    fm.get_jk = _jk

    def _occ(mo_energy, mo_coeff):
        # the filling follows the seed by the per-slot MOM assignment
        # on the CLEAN spin-summed Fock's eigenbasis -- NOT the
        # kernel's columns: the Roothaan effective Fock's low spectrum
        # is MIRROR-DEGENERATE (a known ROHF property of the
        # block-coupled structure; the pairs split at the 1e-7 scale
        # and flip gauge every cycle -- measured: identical Fock, the
        # picks flip from exact to garbage), and the kernel's columns
        # inherit the chaotic in-pair gauge.  The spin-summed Fock
        # (focka+fockb)/2 = h + J - K/2 has isolated physical gaps,
        # so its eigbasis is stable and the anchors' self-overlaps
        # hold.  The picks then MAP to the kernel's columns by max
        # overlap (the assignment owns the configuration; the kernel
        # gauge only hosts it).  NO block-index pools: the subspace
        # Fock's eigvectors interleave freely across the family/window
        # layout (the twins' 4s ranks 9th globally).  The character
        # gates only QUALIFY.  The seedless ion leg fills by pure
        # aufbau within the gates (confined_core_orbitals's own
        # discipline).
        n_mo = mo_coeff.shape[1] if mo_coeff.ndim == 2 else len(mo_energy)
        occ = np.zeros(n_mo)
        w_cl = _w_mask(mo_coeff, p_closed)
        w_op = _w_mask(mo_coeff, p_open)
        cl_allowed = [i for i in range(n_mo) if w_cl[i] > 0.5]
        if has_seed:
            F_can = hc_sub + (vhf_cell[0][0] + vhf_cell[0][1]) * 0.5
            ev_can, U_can = np.linalg.eigh(F_can)
            cl_clean = _mom_assign(prev_cl[0], U_can, list(range(n_mo)),
                                   n_cl_tot, ev_can)
            op_clean = _mom_assign(prev_op[0], U_can,
                                   [i for i in range(n_mo)
                                    if i not in set(cl_clean)],
                                   n_open, ev_can)

            def _map(clean, taken):
                mapped, used = [], set(taken)
                for j in clean:
                    ovl = np.abs(U_can[:, j].conj() @ mo_coeff)
                    free = [i for i in range(n_mo) if i not in used]
                    best = max(free, key=lambda i: ovl[i])
                    mapped.append(best)
                    used.add(best)
                return mapped, used
            cl, used_cl = _map(cl_clean, set())
            occ[np.asarray(cl, dtype=int)] = 2.0
            rest = [i for i in range(n_mo)
                    if occ[i] == 0.0 and 2 * w_op[i] + w_cl[i] > 1.0]
            vop, _ = _map(op_clean, used_cl)
            vop = [i for i in vop if i in set(rest)]
            if len(vop) < n_open:
                free = [i for i in rest if i not in set(vop)]
                vop += sorted(free,
                              key=lambda i: round(float(mo_energy[i]), 6)
                              )[:n_open - len(vop)]
            vop = vop[:n_open]
            occ[np.asarray(vop, dtype=int)] = 1.0
        else:
            cl = sorted(cl_allowed,
                        key=lambda i: round(float(mo_energy[i]), 6))[:n_cl_tot]
            occ[np.asarray(cl, dtype=int)] = 2.0
            rest = [i for i in range(n_mo)
                    if occ[i] == 0.0 and 2 * w_op[i] + w_cl[i] > 1.0]
            vop = sorted(rest,
                         key=lambda i: round(float(mo_energy[i]), 6))[:n_open]
            occ[np.asarray(vop, dtype=int)] = 1.0
        prev_cl[0] = mo_coeff[np.asarray(cl, dtype=int)].T.copy()
        prev_op[0] = mo_coeff[np.asarray(vop[:n_open], dtype=int)].T.copy()
        return occ

    fm.get_occ = _occ

    # the Newton path's gradient tolerance is the certificate's own:
    # the seeds ARE converged solutions (their |g| sits at the free
    # runs' conv_tol_grad ~1e-6, the measured Sc head: 1.08e-6) -- the
    # default sqrt(conv_tol) = 1e-6 threshold misses them by a hair,
    # the needless first step interacts with the mirror-gauge fills
    # and the walk diverges (measured: d2s1 walked from E = the exact
    # anchor to -758.4).  1e-4 recognizes the seed at macro 0; the
    # fresh ion solve converges well inside it.
    if "newton" in conv_accelerator:
        fm.conv_tol_grad = 1e-4

    # the family's ONE acceleration dispatch: plain DIIS on the twins'
    # whole-space landscape is the certified oscillation class (the
    # Stage-2 recipe applies here verbatim -- DIIS holds the seed's
    # neighborhood, the Newton switch closes the tail)
    from atomterm.reference import _accelerated
    from types import SimpleNamespace as _NS
    fm = _accelerated(fm, _NS(conv_accelerator=conv_accelerator), dm0=dm0)
    n_fill = int(round(float(np.sum(fm.mo_occ))))
    c_ao = Q @ fm.mo_coeff
    occ = np.asarray(fm.mo_occ)
    occ_idx = np.flatnonzero(occ > 0)
    splits = {}
    for name, P in (("closed", p_closed), ("open", p_open)):
        w = np.einsum("ai,ab,bi->i", c_ao[:, occ_idx].conj(), P,
                      c_ao[:, occ_idx]).real
        splits[name] = [round(float(x), 4) for x in w]
    report = {
        "n_fam": n_fam, "n_win": n_win, "n_pairs": n_pairs,
        "converged": bool(fm.converged), "e_tot": float(fm.e_tot),
        "n_fill": n_fill, "n_want": int(round(mol.nelectron)),
        "capture": (None if capture is None
                    else [round(float(x), 4) for x in capture]),
        "window_splits": splits,
    }
    if not fm.converged:
        raise RuntimeError(
            f"the subspace SCF did not converge (E = {fm.e_tot:.8f})")
    if n_fill != mol.nelectron:
        raise RuntimeError(
            f"the converged filling holds {n_fill} electrons, want "
            f"{mol.nelectron} -- a lost configuration, not a result")
    return {"e_tot": float(fm.e_tot), "converged": True, "c_ao": c_ao,
            "mo_occ": occ, "report": report}
