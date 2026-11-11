"""Core construction: the purified, exactly rotation-closed core span
and the closed-shell core energy.

PURPOSE
-------
The pur substrate: an exactly rotation-closed replacement for the
canonical core span. Per angular momentum l, the core block is
rebuilt in block form as

    C_pur^{(l)} = E_l ⊗ 1_M

-- one shared radial set E_l (n_r columns) replicated on all M = 2l+1
members. The closure theorem (E ⊗ members is L-invariant for ANY
radial subspace E; README · Methods) then applies to the core itself:
every window orthogonalized against the core -- the (1 - P_core)
projection -- inherits EXACT closure, and the core field is
scalar by construction (its density is exactly m-symmetric).

Why replace the canonical core at all: its MOs carry foreign-l tails
and m-asymmetric same-l content, both of which break the core span's
L-closure, and every window orthogonalized against that span inherits
the defect -- the input-window baseline every CASSCF start carries,
growing with core-window adjacency (the measured classes: the test
pins; the mechanism:
:func:`atomterm.windows.build_full_mo` Notes).
"""

from __future__ import annotations

import numpy as np

from atomterm.windows import (lm_sectors, lowdin, occupied_core_idx,
                              shell_weights)

__all__ = ["purified_core_orbitals", "confined_core_orbitals",
           "closed_shell_core_energy"]

# the closed-core Roothaan's cycle budget (a core span that cannot
# settle is a declaration error, not a patience question)
_CORE_MAXCYCLES = 100


# the core span ------------------------------------------------------

def purified_core_orbitals(mol, mf, ncore: int):
    """The canonical core MOs replaced by an exactly rotation-closed
    core span (the pur substrate).

    The canonical core MOs of an open-shell atom carry foreign-l tails
    and m-ASYMMETRIC same-l content; both break the L-closure of the
    core SPAN, and every window orthogonalized against that span
    inherits the defect -- the input-window baseline, growing with
    core-window adjacency (measured classes: the test pins).

    Notation, per angular momentum l: the shell has M = 2l+1 member
    sectors, each an identical radial chain of n_r functions with
    overlap S_r; n_l canonical core MOs are dominated by l. Build

        C^m = C[rows of member m, the n_l MOs]   (n_r x n_l, per m)
        R   = [C^1  C^2  ...  C^M]               (n_r x M·n_l)
        G   = R^† S_r R                          (M·n_l x M·n_l)

    In the m-symmetric ideal (every member carrying the same radials,
    C^m = C0 for all m) the Gram collapses to

        G = Σ_m C^m† S_r C^m = M · C0^† S_r C0,

    so each of the k_l = n_l/M occupied radials contributes ONE
    dominant eigenvalue ≈ M · (its S-weight), and everything beyond
    k_l is tail. The construction reads exactly that spectrum:

        Λ_k, Q_k = the top-k_l eigenpairs of G     (the beyond-k
        spectrum is asserted tail-level -- the tail level is both the
        selection and the certificate)
        E = R Q_k Λ_k^{-1/2}                       (n_r x k_l,
        S_r-orthonormal: E^† S_r E = I)
        C_pur[member m of l] = E  for EVERY m      (k_l x M columns,
        m-complete by construction)

    One global Lowdin then orthonormalizes C_pur (S = 0 between
    different l keeps l-purity; within-l mixing is span-preserving).

    Parameters
    ----------
    mol : the built Mole (the shell walker's substrate).
    mf : converged reference (its ``mo_coeff`` and ``get_ovlp``).
    ncore : number of core orbitals (the closed block;
        :func:`atomterm.windows.occupied_core_idx` selects it).

    Returns
    -------
    (C, contamination) : ``(nao, ncore)`` S-orthonormal columns, and
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
    TRAPS: solve in the COLUMN space (``G = R^† S_r R``); the
    row-space S-generalized eigenproblem ``R R^† v = λ S_r v``
    returns S^{-1}-warped radials. And take R from the FULL l-rows of
    C -- canonical core MOs are NOT shell-confined (C 2s: 66% of
    S-weight off its shell), so per-shell projection loses the tail.
    """
    s = mf.get_ovlp()
    c_can = np.asarray(
        mf.mo_coeff[:, occupied_core_idx(mf, ncore)])
    sectors = lm_sectors(mol)

    # dominant l per core MO (lumped-l weights partition unity exactly;
    # the foreign-l weight is the tiny numerical tail)
    l_names, w_l = shell_weights(mol, c_can, s)
    dom_l = [l_names[i] for i in w_l.argmax(axis=1)]
    for i, wmax in enumerate(w_l.max(axis=1)):
        assert wmax > 0.95, (i, wmax)  # a mixed-l core MO = wrong core

    c_pur = np.zeros((mol.nao_nr(), ncore), dtype=c_can.dtype)
    contamination = np.zeros(ncore)
    col = 0
    for k_col, l_name in enumerate(l_names):
        l_val = "spdfgh".index(l_name)
        mos = [i for i, d in enumerate(dom_l) if d == l_name]
        if not mos:
            continue
        members = sorted(m for (ll, m) in sectors if ll == l_val)
        k = len(mos) // len(members)   # k_l = n_l / M (the m-set count)
        assert k * len(members) == len(mos), (
            f"core MOs of l={l_name} do not form complete m-sets: {len(mos)}")
        chain = sectors[(l_val, members[0])]
        s_rad = s[np.ix_(chain, chain)]          # S_r (identical per member)

        # foreign-l contamination diagnostic
        for i in mos:
            contamination[i] = 1.0 - w_l[i, k_col]

        # R = [C^1 ... C^M]: the member-resolved radial blocks (n_r x M·n_l)
        rr = np.hstack([c_can[sectors[(l_val, m)], :][:, mos]
                        for m in members])
        # G = R^† S_r R: top-k eigenvalues = the M-scaled radials,
        # beyond-k = tail
        ev, q = np.linalg.eigh(rr.conj().T @ s_rad @ rr)
        order = np.argsort(ev)[::-1]
        top = order[:k]
        assert ev[top[-1]] > 1e-8, (l_name, ev)
        if len(ev) > k:  # discarded directions are tails, not radials
            assert ev[order[k]] < 1e-2 * ev[top[-1]], (l_name, ev)
        v_all = rr @ q[:, top] @ np.diag(ev[top] ** -0.5)  # E = R Q_k Λ_k^{-1/2}
        for j in range(k):
            for m in members:  # C_pur^m = E for EVERY member (m-complete)
                c_pur[sectors[(l_val, m)], col] = v_all[:, j]
                col += 1
    assert col == ncore

    # one global Loewdin (S = 0 between different l -> l-purity survives;
    # within-l mixing is span-preserving)
    return lowdin(c_pur, s), contamination


# the closed core on a core span --------------------------------------

def confined_core_orbitals(mf, c_span: np.ndarray,
                           n_pairs: int) -> tuple[np.ndarray, dict]:
    """The closed core ON a core span: the span's OWN best closed-shell
    pairs -- a Loewdin-span-basis closed-shell Roothaan, self-consistent
    under the nuclear field and the core's own two-body field.

    A shell-span core declares its CHEMISTRY by the electron count
    (``n_pairs``) and its AO realization by the span; on a
    multi-column realization (n_pairs pairs over M > n_pairs columns)
    the frozen core must be selected INSIDE the span, and the
    span-truncated projection of the SCF core is NOT that selection:
    the S-weight capture is a NORM statement blind to the nuclear
    cusp (a tiny-norm, Eh-scale region -- measured on Sc-unc: 99.99%
    capture with +13 Eh core-energy error). The span's own minimum
    is: diagonalize, fill the lowest ``n_pairs`` doubly, rebuild the
    closed-shell Fock, iterate -- converging in a few cycles to the
    closed-shell (ion-like) core the span can host, monotone in the
    span, and reading NO converged orbitals (the core stays
    draw-free). The span's MEMBERSHIP carries the shell grammar: the
    lowest pairs of the core's family span are its declared shells,
    so a core that skips a shell must leave that shell's bands out of
    its realization.

    The degenerate shortcut ``n_pairs == M`` returns the span itself:
    every column then IS one doubly-occupied orbital (the certified
    1:1 raw-shell core -- this path reproduces it bit-exactly).

    Parameters
    ----------
    mf : the reference (``get_hcore``; ``mol`` for the eri). No
        converged orbitals are read.
    c_span : (nao, M) S-orthonormal span columns (the core's raw-AO
        span, :func:`atomterm.cells._lane_core`).
    n_pairs : the core's doubly-occupied pair count (the CHEMISTRY:
        :func:`atomterm.windows.closed_electrons` over the core's
        shells, halved).

    Returns
    -------
    (c_core, report) : ``(nao, n_pairs)`` S-orthonormal columns, and
        the mini-SCF's convergence report ``{"cycles", "comm"}``
        (``comm`` = the final density-Fock commutator norm -- the
        self-consistency certificate; ``{"cycles": 0, "comm": 0.0}``
        on the 1:1 path).

    Raises
    ------
    ValueError
        If ``n_pairs`` exceeds the span width, or the Roothaan loop
        does not converge (``max_cycles``) -- a span the closed-shell
        problem cannot settle in is a declaration error.

    Notes
    -----
    The converged columns are occupied eigenvectors of their own
    Fock, so the core's bookkeeping :func:`closed_shell_core_energy`
    over them is the span's closed-shell variational value; on the
    whole AO space the construction reproduces the (2·n_pairs)-ion
    RHF exactly (pinned in the tests).
    """
    from pyscf import ao2mo

    m = c_span.shape[1]
    if n_pairs > m:
        raise ValueError(
            f"{n_pairs} closed pairs over a {m}-column span -- the core "
            "declares more electrons than its realization holds columns")
    if n_pairs == m:
        return c_span, {"cycles": 0, "comm": 0.0}

    h = c_span.conj().T @ mf.get_hcore() @ c_span
    eri = ao2mo.restore(1, ao2mo.kernel(mf.mol, c_span), m)
    f = h.copy()
    comm = float("inf")
    occ = None
    cycles = 0
    for cycles in range(1, _CORE_MAXCYCLES + 1):
        _, u = np.linalg.eigh(f)
        occ = u[:, :n_pairs]
        d = 2.0 * (occ @ occ.conj().T)          # the closed-shell density
        f = (h
             + np.einsum("rs,pqrs->pq", d, eri, optimize=True)
             - 0.5 * np.einsum("rs,prqs->pq", d, eri, optimize=True))
        # the stationarity measure: the total energy settles long
        # before the orbitals do (quadratic vs linear), so the break
        # reads the DENSITY-Fock commutator itself
        comm = float(np.abs(f @ d - d @ f).max())
        if comm < 1e-10:
            break
    else:
        raise ValueError(
            f"the closed-core Roothaan did not settle in {_CORE_MAXCYCLES}"
            f" cycles (last |[F, D]| = {comm:.2e}) -- the core span "
            "declaration needs revisiting")
    return c_span @ occ, {"cycles": cycles, "comm": comm}


# the core energy ----------------------------------------------------

def closed_shell_core_energy(mf, c_core: np.ndarray) -> float:
    """Energy of a doubly-occupied core span in the raw mf field.

    The pur-core bookkeeping of the campaigns (promoted from the scan
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
        bookkeeping constant of the core, not a variational energy.
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
