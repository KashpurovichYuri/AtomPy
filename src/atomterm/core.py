"""Core construction: the purified, exactly isotropic core span and the
closed-shell core energy.

The pur substrate (fact 20) rebuilds the closed core span as
E_l \u2291 {all members of l} per angular momentum -- \u00a72.3's closure
theorem applied to the core itself -- so that the start window's
(1 - P_core) projection inherits EXACT closure, and the fold field is
scalar by construction.
"""

from __future__ import annotations

import numpy as np

from atomterm.windows import (lm_sectors, lowdin, occupied_core_idx,
                              shell_weights)

__all__ = ["purified_s_core", "purified_core_orbitals",
           "closed_shell_core_energy"]


def purified_s_core(mol, mf, n_core: int):
    """Canonical core MOs projected onto the s-AO subspace (l = 0 exactly).

    The converged "s" MOs of an open-shell atom carry a tiny p/d tail
    (measured on C/cc-pVDZ: 1s 1.9e-7, 2s 6.3e-6 weight deficit) which
    breaks two preconditions at once -- orthonormality (<MO|2p_AO>
    ~ 1e-3, so second quantization is ill-defined on such columns) and
    rotational closure ([H, Lz] fails at ~1e-4 instead of machine
    precision). Projection + a Loewdin inside the s subspace repairs
    both at second-order cost (~1e-10 Eh) -- the orbital analog of the
    density-space spherical average. Full causal chain: fact 9.

    Returns ``(C, contamination)``: (nao, n_core) orthonormal columns
    and the dropped fraction ``1 - <c_s|S|c>`` per canonical MO (the SCF
    symmetry-breaking diagnostic; MO weights live in the AO metric,
    c^dag S c = 1, not the Euclidean norm).
    """
    import re

    labels = [line.split()[-1] for line in mol.ao_labels()]
    s_idx = [i for i, lbl in enumerate(labels) if re.fullmatch(r"\d+s", lbl)]

    c_can = mf.mo_coeff[:, occupied_core_idx(mf, n_core)]
    c_s = c_can[s_idx, :]                                # drop the p/d tail
    # MO normalization lives in the AO metric (c^dag S c = 1), so the
    # s-subspace weight is c_s^dag S_ss c_s -- not the Euclidean norm.
    s_ov = mf.get_ovlp()[np.ix_(s_idx, s_idx)]
    weight = np.einsum("ij,ij->j", c_s, s_ov @ c_s)
    contamination = np.clip(1.0 - weight, 0.0, None)     # fp noise can dip < 0
    c_s = c_s / np.sqrt(weight)                          # renormalize the rest

    c_s = lowdin(c_s, s_ov)

    c = np.zeros((mol.nao_nr(), n_core), dtype=c_s.dtype)
    c[s_idx, :] = c_s
    return c, contamination


def purified_core_orbitals(mol, mf, n_core: int):
    """Canonical core MOs replaced by an exactly rotation-closed core span.

    Generalizes :func:`purified_s_core` to p/d cores (Si [Ne]+3s, Sc
    [Ar]). The canonical core MOs of an open-shell atom carry foreign-l
    tails and m-ASYMMETRIC same-l content; both break the L-closure of
    the core SPAN, and every window orthogonalized against that span
    inherits the defect (the fact-17 input-window baseline, ~1e-3 Si /
    ~1e-2 Sc). Closure theorem (E (x) members is L-invariant for any
    radial subspace E): fact 20, README §5.3. Construction per l:

    1. member-resolved radial blocks ``R = [r_{i,m}]`` (each column =
       MO i's coefficients on member m's radial chain);
    2. E_l = top-k_l span of the S-metric column Gram ``R^T S_rad R``,
       S-orthonormalized -- eigenvalues = S-weight x member multiplicity
       (~2l+1 per occupied radial), so the beyond-k spectrum is the tail
       content, asserted as such;
    3. every radial of E_l replicated on ALL member sectors of l
       (k_l x (2l+1) columns -- m-complete by construction).

    One global Loewdin then orthonormalizes the set (S = 0 between
    different l keeps l-purity; within-l mixing is span-preserving).
    TRAPS (fact 20): the S-generalized eigenproblem of ``R R^T`` returns
    S^{-1}-warped radials -- the column Gram is the correct object; and
    canonical core MOs are NOT shell-confined (C 2s: 66% of S-weight off
    its shell) -- per-shell projection is wrong. For l = 0 the
    construction reduces exactly to :func:`purified_s_core`.

    Returns ``(C, contamination)``: (nao, n_core) S-orthonormal columns
    and per-MO ``1 - own-l weight`` (the foreign-l SCF tail diagnostic).
    """
    s = mf.get_ovlp()
    c_can = np.asarray(
        mf.mo_coeff[:, occupied_core_idx(mf, n_core)])
    sectors = lm_sectors(mol)

    # dominant l per core MO (lumped-l weights partition unity exactly;
    # the foreign-l weight is the tiny numerical tail, fact 9)
    l_names, w_l = shell_weights(mol, c_can, s)
    dom_l = [l_names[i] for i in w_l.argmax(axis=1)]
    for i, wmax in enumerate(w_l.max(axis=1)):
        assert wmax > 0.95, (i, wmax)  # a mixed-l core MO = wrong core

    c_pur = np.zeros((mol.nao_nr(), n_core), dtype=c_can.dtype)
    contamination = np.zeros(n_core)
    col = 0
    for k_col, l_name in enumerate(l_names):
        l_val = "spdfgh".index(l_name)
        mos = [i for i, d in enumerate(dom_l) if d == l_name]
        if not mos:
            continue
        members = sorted(m for (ll, m) in sectors if ll == l_val)
        k = len(mos) // len(members)   # radial copies (n_mos / m-multiplicity)
        assert k * len(members) == len(mos), (
            f"core MOs of l={l_name} do not form complete m-sets: {len(mos)}")
        chain = sectors[(l_val, members[0])]
        s_rad = s[np.ix_(chain, chain)]

        # foreign-l contamination diagnostic
        for i in mos:
            contamination[i] = 1.0 - w_l[i, k_col]

        # dominant radial subspace E_l: top-k of the S-metric column Gram
        # of the member-resolved radial blocks
        rr = np.hstack([c_can[sectors[(l_val, m)], :][:, mos]
                        for m in members])
        ev, q = np.linalg.eigh(rr.conj().T @ s_rad @ rr)
        order = np.argsort(ev)[::-1]
        top = order[:k]
        assert ev[top[-1]] > 1e-8, (l_name, ev)
        if len(ev) > k:  # discarded directions are tails, not radials
            assert ev[order[k]] < 1e-2 * ev[top[-1]], (l_name, ev)
        v_all = rr @ q[:, top] @ np.diag(ev[top] ** -0.5)  # S-orthonormal
        for j in range(k):
            for m in members:  # shared radial on EVERY member (m-complete)
                c_pur[sectors[(l_val, m)], col] = v_all[:, j]
                col += 1
    assert col == n_core

    # one global Loewdin (S = 0 between different l -> l-purity survives;
    # within-l mixing is span-preserving)
    return lowdin(c_pur, s), contamination



def closed_shell_core_energy(mf, c_core: np.ndarray) -> float:
    """E of a doubly-occupied core in the raw mf field (the pur-fold
    bookkeeping of the campaigns; promoted from the scan scripts)."""
    from pyscf import ao2mo

    hc = c_core.T @ mf.get_hcore() @ c_core
    n = c_core.shape[1]
    eri = ao2mo.restore(1, ao2mo.kernel(mf.mol, c_core), n)
    e = 2.0 * float(np.trace(hc))
    for i in range(n):
        e += float(np.sum(2.0 * np.diag(eri[i, i])
                          - eri[i, :, :, i].diagonal()))
    return e
