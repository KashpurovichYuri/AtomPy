"""Core construction: the purified, exactly rotation-closed core span
and the closed-shell core energy.

PURPOSE
-------
The pur substrate (fact 20): rebuild the closed core span as
E_l ⊗ {all members of l} per angular momentum -- the closure theorem
(E ⊗ members is L-invariant for ANY radial subspace E; README ·
Methods) applied to the core itself -- so that every start window's
(1 - P_core) projection inherits EXACT closure and the fold field is
scalar by construction. This is the fix for the fact-17 input-window
baseline: canonical core MOs carry foreign-l tails and m-asymmetric
same-l content, and any window orthogonalized against them inherits
the defect.
"""

from __future__ import annotations

import numpy as np

from atomterm.windows import (lm_sectors, lowdin, occupied_core_idx,
                              shell_weights)

__all__ = ["purified_core_orbitals", "closed_shell_core_energy"]


def purified_core_orbitals(mol, mf, n_core: int):
    """The canonical core MOs replaced by an exactly rotation-closed
    core span (the pur substrate; fact 20).

    The canonical core MOs of an open-shell atom carry foreign-l tails
    and m-ASYMMETRIC same-l content; both break the L-closure of the
    core SPAN, and every window orthogonalized against that span
    inherits the defect (the fact-17 input-window baseline, ~1e-3 Si /
    ~1e-2 Sc). The construction, per angular momentum l:

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

    Parameters
    ----------
    mol : the built Mole (the shell walker's substrate).
    mf : converged reference (its ``mo_coeff`` and ``get_ovlp``).
    n_core : number of core orbitals (the closed block;
        :func:`atomterm.windows.occupied_core_idx` selects it).

    Returns
    -------
    (C, contamination) : ``(nao, n_core)`` S-orthonormal columns, and
        per-MO ``1 - own-l weight`` of the canonical core it replaces
        (the foreign-l SCF tail diagnostic).

    Raises
    ------
    AssertionError
        On a mixed-l core MO (dominant weight <= 0.95 = wrong core),
        core MOs not forming complete m-sets, or a discarded radial
        direction that is not tail-level -- each means the input core
        violates the construction's contract.

    Notes
    -----
    TRAPS (fact 20): the S-generalized eigenproblem of ``R R^T``
    returns S^{-1}-warped radials -- the column Gram is the correct
    object; and canonical core MOs are NOT shell-confined (C 2s: 66%
    of S-weight off its shell) -- per-shell projection is wrong.
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
    """Energy of a doubly-occupied core span in the raw mf field.

    The pur-fold bookkeeping of the campaigns (promoted from the scan
    scripts): with ``h = C^dag h_core C`` and the chemist eri over the
    core columns,

        E = 2 Σ_i h_ii + Σ_ij (2 (ii|jj) − (ij|ji))
          = 2 tr h + Σ_ij (2 J_ij − K_ij).

    Parameters
    ----------
    mf : the reference object (``get_hcore``; ``mol`` for the eri).
    c_core : (nao, n) S-orthonormal core columns.

    Returns
    -------
    float
        The closed-shell (doubly-occupied) energy of the span -- a
        bookkeeping constant of the fold, not a variational energy.
    """
    from pyscf import ao2mo

    hc = c_core.T @ mf.get_hcore() @ c_core
    n = c_core.shape[1]
    eri = ao2mo.restore(1, ao2mo.kernel(mf.mol, c_core), n)
    e = 2.0 * float(np.trace(hc))
    for i in range(n):
        e += float(np.sum(2.0 * np.diag(eri[i, i])
                          - eri[i, :, :, i].diagonal()))
    return e
