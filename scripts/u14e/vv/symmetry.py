"""The D2h column algebra on a symmetry-tagged mol: the metric
projectors, the irrep cleanup, the per-irrep orthonormalization, the
core/irrep intersection."""

import numpy as np  # noqa: E402

def _irrep_projectors(mol, S):
    """``[(D_k, G_k^{-1} D_k^T S)]`` per irrep block -- the metric
    projector pieces, ``P_k v = D_k @ (Ginv_dTS @ v)``.  The irrep
    spans are mutually S-orthogonal and partition the AO space, so
    ``sum_k ||P_k v||^2_S = ||v||^2_S`` whatever the blocks' internal
    gauge (pyscf keeps raw AO selectors in single-irrep blocks,
    max|D^T S D - I| = 0.25 on the sto-3g probe)."""
    projs = []
    for d in mol.symm_orb:
        d = np.asarray(d)
        G = d.T @ S @ d
        projs.append((d, np.linalg.solve(G, d.T @ S)))
    return projs


def symmetrize_columns(mol, S, c, tol: float = 1e-4):
    """Re-project each column onto its DOMINANT D2h irrep -- the
    disclosed CLEANUP for the tiny irrep tails a core projection's
    conditioning leaves (``orb_irreps``/``label_orb_symm`` refuses
    past ~1e-6); the core projection itself is physics and is never
    undone here.  The projection is the METRIC projector
    ``P_k = D_k (D_k^T S D_k)^{-1} D_k^T S`` per irrep block: the
    irrep spans are mutually S-orthogonal and partition the AO space,
    so ``sum_k ||P_k v||^2_S = ||v||^2_S`` exactly and the impurity
    ``1 - max_k ||P_k v||^2_S / ||v||^2_S`` sits in ``[0, 1]``
    whatever the blocks' internal gauge (a coordinate-only
    projection overcounts when ``D_k`` is not S-orthonormal -- the
    defect the certification caught).  The impurity is measured on
    the S-normalized column, so it is scale-free; already-pure
    columns project to themselves (the identity on their span) and
    come back just normalized.  Returns ``(c_fixed, impurities)`` --
    the impurities are the per-column BEFORE values, the caller
    discloses them; raises past ``tol`` (a real mixing class, not
    the projection-conditioning dirt this cleanup owns)."""
    projs = _irrep_projectors(mol, S)
    fixed = np.empty_like(c)
    imp = []
    for j in range(c.shape[1]):
        v = np.asarray(c[:, j], dtype=float)
        n2 = float(v @ S @ v)
        best_w, best_n = None, -1.0
        for d, Ginv_dT in projs:
            w = d @ (Ginv_dT @ v)
            n_k = float(w @ S @ w)
            if n_k > best_n:
                best_w, best_n = w, n_k
        imp.append(1.0 - best_n / n2)
        fixed[:, j] = best_w / np.sqrt(best_n)
    worst = max(imp)
    if worst > tol:
        raise ValueError(
            f"a column carries irrep impurity {worst:.2e} past the "
            f"tolerance {tol:.1e} -- a real mixing class, not the "
            f"projection-conditioning dirt this cleanup owns")
    return fixed, imp


def dominant_irreps(mol, S, c):
    """Per column, the id of the irrep span holding the most of it
    (the metric-projector weights); :func:`symmetrize_columns`
    re-projects onto these."""
    projs = _irrep_projectors(mol, S)
    ids = []
    for j in range(c.shape[1]):
        v = np.asarray(c[:, j], dtype=float)
        best_k, best_n = 0, -1.0
        for k, (d, Ginv_dTS) in enumerate(projs):
            w = d @ (Ginv_dTS @ v)
            nk = float(w @ S @ w)
            if nk > best_n:
                best_k, best_n = k, nk
        ids.append(best_k)
    return ids


def orthonormal_per_irrep(mol, S, act):
    """S-orthonormalize the act block WITHIN each dominant-irrep
    group -- the pre-assembly treatment for blocks whose columns are
    irrep-pure but not mutually S-orthogonal (the MOs-window class:
    the replaced s column overlaps the f/d columns at ~1e-2, and
    ``build_full_mo``'s global Loewdin would mix the pure columns
    ACROSS irrep blocks -- the assembled block then fails pyscf's
    orb labeling at the ~1e-04 scale).  Per-group Loewdin never
    leaves a group, so every column's purity is preserved exactly;
    the groups are mutually S-orthogonal, so the block comes out
    S-orthonormal overall and the assembly's global Loewdin becomes
    the identity.  Returns (act_fixed, worst_purity_loss) -- the loss
    is the largest per-column dominant-block weight defect (the
    caller discloses it; ~1e-16 when the input is already clean)."""
    projs = _irrep_projectors(mol, S)
    ids = dominant_irreps(mol, S, act)
    fixed = np.empty_like(act)
    worst = 0.0
    for k in sorted(set(ids)):
        cols = [j for j, kid in enumerate(ids) if kid == k]
        sub = np.asarray(act)[:, cols]
        g = sub.T @ S @ sub
        gv, U = np.linalg.eigh(g)
        sub = sub @ U @ np.diag(gv ** -0.5) @ U.T
        fixed[:, cols] = sub
        d, Ginv_dTS = projs[k]
        for q in range(sub.shape[1]):
            w = d @ (Ginv_dTS @ sub[:, q])
            n2 = float(sub[:, q] @ S @ sub[:, q])
            worst = max(worst, 1.0 - float(w @ S @ w) / n2)
    return fixed, worst


def purify_window(mol, S, c_core, c, tol: float = 1e-10):
    """Converge a window block to the INTERSECTION of the frozen-core
    S-orthogonalization and the dominant-irrep purity, in closed form.

    The naive iterate (project core, project irrep, renormalize)
    converges only geometrically on a core-lean column (measured at
    (4,15): rate ~0.36 at 92% core-lean -- stalled at 5e-04 after 8
    passes; a single pass leaves a 2.7e-07 tail past pyscf's orb
    labeling).  The closed form removes the core's content
    within each dominant-irrep block exactly: the block subtracts
    ``B_k = P_k(c_core)`` -- ALL core columns projected into
    ``span(D_k)`` -- so no assignment of core columns to blocks is
    needed (the lane core's global orthonormalization mixes across
    irrep blocks; a per-assigned-block form stalls at 5.3e-04).

    Returns ``(c_fixed, report)`` with ``report = (passes,
    core_res, irrep_imp, ortho_dev)``; raises on genuinely-mixed
    input (via :func:`symmetrize_columns`' contract), on a vanished
    column (the column lies INSIDE the frozen core's span -- it has
    no core-orthogonal remainder at all; the measured (4,15) class:
    the C1 reference's m-mixed occupied core, S-projected onto the
    declared shell span, contributes the promoted band's own slots),
    or when the verification residues miss ``tol``.  THE CONTRACT:
    the caller passes the ASSEMBLED block (the pipeline's form --
    ``build_full_mo``'s output, already core-orthogonal); a
    genuinely-core-lean RAW column is refused -- core-projecting it
    is the assembly's job, and a global pre-projection here would
    SYNTHESIZE irrep content against a cross-irrep-mixed core (the
    defect the certification caught: the remainder picks up the
    core's other-block directions at the correlation scale).  THE
    CORE NEED NOT BE IRREP-PURE: the lane core's global
    orthonormalization mixes across irrep blocks (measured:
    per-assigned-block subtraction stalls at a 5.3e-04 core
    residue), so each block removes the core's COMPLETE content
    within its span, ``B_k = P_k(c_core)`` -- assignments
    irrelevant."""
    v, imp = symmetrize_columns(mol, S, c)
    ids = dominant_irreps(mol, S, v)
    projs = _irrep_projectors(mol, S)
    fixed = np.empty_like(v)
    for k in sorted(set(ids)):
        cols = [j for j, kid in enumerate(ids) if kid == k]
        d, Ginv_dTS = projs[k]
        sub = v[:, cols]
        B = d @ (Ginv_dTS @ c_core)     # ALL core columns, in span(D_k)
        Gb = B.T @ S @ B
        gb, Ub = np.linalg.eigh(Gb)
        keep = gb > 1e-10 * gb.max()
        Bk = B @ Ub[:, keep] / np.sqrt(gb[keep])
        sub = sub - Bk @ (Bk.T @ S @ sub)   # the block's complete core
        sub = d @ (Ginv_dTS @ sub)          # the irrep cleanup: no-op
        norms = np.sqrt(np.diag(sub.T @ S @ sub))
        if norms.min() < 1e-8:
            raise VanishedColumnError(
                f"a window column vanished into the core (post-"
                f"subtraction S-norm {norms.min():.1e}) -- the column "
                f"lies inside the frozen core's span and has no core-"
                f"orthogonal remainder")
        sub = sub / norms
        g = sub.T @ S @ sub
        gv, U = np.linalg.eigh(g)
        if gv.min() <= 1e-10 * max(gv.max(), 1.0):
            raise ValueError(
                f"the window block is rank-deficient in irrep {k} "
                f"(smallest Gram eigenvalue {gv.min():.1e}) -- "
                f"dependent columns, not a window")
        sub = sub @ U @ np.diag(gv ** -0.5) @ U.T
        for q, j in enumerate(cols):
            fixed[:, j] = sub[:, q]
    core_res = float(np.linalg.norm(c_core.T @ S @ fixed, axis=0).max())
    _, imp_v = symmetrize_columns(mol, S, fixed, tol=1.0)
    ortho_dev = float(np.max(np.abs(
        fixed.T @ S @ fixed - np.eye(fixed.shape[1]))))
    if core_res > tol or max(imp_v) > tol or ortho_dev > 1e-8:
        raise ValueError(
            f"the closed-form intersection missed its tolerances "
            f"(core residue {core_res:.1e}, irrep impurity "
            f"{max(imp_v):.1e}, orthonormality {ortho_dev:.1e})")
    return fixed, {"passes": 1, "core_res": core_res,
                   "irrep_imp": float(max(imp_v)), "ortho_dev": ortho_dev}


class VanishedColumnError(ValueError):
    """A window column vanished into the core inside
    :func:`purify_window` -- the column lies inside the frozen core's
    span (no core-orthogonal remainder exists for it; the measured
    (4,15) class: the C1 reference's m-mixed occupied core, projected
    onto the declared shell span, spans the promoted band's own
    slots)."""
