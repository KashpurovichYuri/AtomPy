# Ticket 13: generalized core purification (the closed core span)

Status: RESOLVED (2026-10-15, fact 20; `purified_core_orbitals` +
the `*_pur` ladder rungs + the sa_decomposition purified pass).

## Question

`purified_s_core` handles s-type cores. How to purify p-core (Si 2p⁶) or
future d-core (Cu 3d¹⁰) orbitals -- and can the purified core span make
the CASSCF input window exactly rotation-closed?

## Resolution (as implemented -- the spec's premise was corrected mid-air)

The approved spec said "project each core MO onto its OWN (n,l) shell and
the closed core spans each FULL AO subspace". GATE (A) REFUTED THE
PREMISE within minutes: canonical core MOs are NOT own-n confined -- C's
2s MO carries 66% of its S-weight on other radials (same-l radial mixing
is physical), so own-shell truncation destroys the reference (density
change 0.77). The closure theorem that actually holds is sharper and
simpler:

  For a one-center basis the l-chain decomposes into MEMBER sectors with
  identical radial chains (S = I_m (x) S_rad exactly), and L acts ONLY on
  the member index. Hence E (x) members is L-invariant for ANY radial
  subspace E: shared radial per closed shell is ALL that closure needs;
  full AO shell subspaces are unnecessary (and physically wrong).

Construction (`purified_core_orbitals(mol, mf, n_core)`):
1. dominant l per core MO (lumped-l S-weights partition unity; asserted
   > 0.95);
2. per l, the member-resolved radial blocks R = [r_{i,m}] of those MOs;
3. E_l = top-k span of the S-metric COLUMN Gram R^T S_rad R (eigenvalues
   count S-weight x member multiplicity ~ 2l+1 per occupied radial --
   the discarded spectrum is asserted tail-level: selection AND
   certificate in one spectrum). METRIC TRAP recorded: the S-generalized
   eigenproblem of R R^T returns S^-1-warped radials (wrong span) --
   the column Gram is the correct object;
4. every radial of E_l replicated on ALL member sectors (m-complete by
   construction), one global Loewdin (S = 0 between different l keeps
   l-purity). k = n_mos(l)/(2l+1) (m-sets asserted complete).
   For l = 0 this reduces exactly to `purified_s_core`.

Verified gates (C / Si / Sc):
(a) reference shift = the discarded core polarization only: spherically
    averaged density change 4.2e-6 / 3.6e-6 / 1.05e-4; energy (same
    projected window, canonical vs purified core) +3.0 uEh / +128 uEh /
    +674 uEh -- C's equals the fact-12 avg-vs-raw shift exactly. The
    purified core IS the orbital realization of the scalar core the fold
    already uses: its density is m-scalar EXACTLY (equal member blocks,
    zero cross-member coherences, anisotropy 1e-16, raw fold == averaged
    fold).
(b) HEADLINE: input windows (build_full_mo AND radial_basis) close at
    machine precision -- Si 3.4e-3/4.5e-3 -> 1.1e-15/1.8e-15, Sc
    5.7e-3/2.1e-2 -> 6.6e-15/8.7e-15. The projector onto E (x) members
    commutes with L, so the window projection leaks nothing.
(c) THE DELIVERABLE -- Sc CASSCF/SA-CASSCF with unbroken Lz and L^2:
    relax_pur 1.3e-13; sa_relax_pur (NEW RUNG: SA + core_virtual_only
    mask -- window frozen, core relaxes under the ensemble) 6.7e-14 /
    3.6e-13, full 90-state free-ion term table, max|dE| 5e-8. Free SA:
    C 0.0 exactly, Si 3.3e-16 (vs 4e-6 healing on canonical) with the
    full fact-18 gains -- the radial direction's restoring force makes
    SA's stationary point the symmetric window. On Sc the member gauge
    is exactly FLAT (fact 18d), so sa_radial_pur stops at a ~1e-6 window
    tilt (comm 5.7e-7/1.8e-5) that no conv_tol tightening removes
    (1e-7 -> 1e-12 leaves [P,Lp] ~ 1e-6, energy identical) -- an
    optimizer floor with no restoring force, not a baseline. The free-sa
    4s<->3p swap attractor SURVIVES the purified substrate (window
    p = 1/6 exactly again): it is an optimizer basin -- the fact-17
    attribution question closed. Single-root casscf_radial_pur keeps the
    fact-16 member-gauge chaos (C: [P,Lp] 1.1 with [H,Lz] 0.0 -- the
    blind-spot demo).
(d) purification energy effect = (a)'s numbers; pur-rung energies carry
    the fact-14 window-definition offset vs raw-fold CASCI rows
    (Sc +88.7 mEh) -- compare within the pur family; relax-family energy
    SIGNS vs the unpolarized CASCI are fold bookkeeping (our averaged
    fold discards the polarization the optimizer bought; max|dE| ~ 1e-3
    is the detector) -- bounded in tests, not pinned.

Phase 4 (post-hoc window repair): NOT NEEDED -- free SA is machine-exact
on p2 atoms on this substrate.

Validation: 5 new tests (83 total, suite green, ~4.5 min):
`test_purified_core_orbitals_exact_closed_span` (orthonormality, one-l
purity per column, lz AND lperp certificates, m-scalar density,
contamination/density bands), `test_purified_substrate_closes_the_input_
windows` (gate (b) on both window constructions, all three atoms),
`test_purified_rungs_exact_closure` (C/Si rung bands + gate (d)),
`test_sc_purified_rungs_unbroken_lz_and_l2` (the deliverable pins,
swap-basin branching). Decomposition artifacts regenerated with the
purified before/after baselines (results/*_sa_decomposition.md).
