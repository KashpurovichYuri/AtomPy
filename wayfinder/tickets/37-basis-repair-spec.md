# Ticket 37 — the basis-repair bench, the @s wave, the subspace-ROHF: the executable spec

**Status**: ACTIVE (the spec finalized 2026-11-12; execution per the Progress table below)
**Depends on**: ticket 36 (the runs wave B→C→D→E→A = Stage 0; approved; independent of this spec's seams)
**Origin**: the 2026-11-12 plan session — the user's "generally agree"; the directives: the to-spec form + progress tracking (this file), **64 OMP threads for the bench**, finalize before proceed.

## Problem Statement

The uncontracted twins (Sc-unc / Ti-unc) carry the measured fresh-SCF disease: every
fresh forced-occupation SCF on these bases diverges (E ~ +1e7 Eh, Newton and DIIS
alike — the d3 census on Sc-unc, every arm; map Issues, 2026-11-11). The disease
blocks: the prior port's `@s` half (the series' non-head configurations —
`cases._unc_s_block` guards the windows as pending rows), the twins' UHF anchor
rows, and any future fresh reference on these bases. Beyond the repair, the
campaign lacks the mean-field certificate that says whether the distilled radial
contractions (`@r1/@s1/@s2`) cost energy at the ROHF level — the certificate whose
adoption decision un-gates the pending optimizer classes (the prior treatments'
re-homing, the twins' primary legs).

## Solution

Two composed repair tracks — P (projection) and C (conditioning) — developed and
certified in a sandbox, then integrated; the twins' `@s` legs then run on whichever
source class won; then the subspace-ROHF exploration (`run_folded`'s mirror)
produces the per-contraction ΔE-vs-anchor certificate, Sc first. Every stage is
gated, tracked below, and loses nothing to compaction (this file + message.txt are
the record).

## User Stories

1. As the campaign, I want the S spectra of both twins' bases measured per
   l-block, so that the conditioning design rests on a measured gap, not an
   assumption — and the design is revisited before any cut if the spectrum decays
   with no gap.
2. As the campaign, I want ONE cross-basis projector primitive, so that the seed
   use and the density-direct use share certified algebra instead of duplicating
   it.
3. As the suite, I want the B=A identity and the B⊇span(A) exactness anchors, so
   that the projector's conventions are pinned with zero parameters and no
   reference calculation.
4. As the suite, I want the near-singular stress with and without the internal
   truncation, so that the conditioning's contribution is demonstrated, not
   asserted.
5. As the physicist, I want the known-answer transfer (the contracted heads'
   occupieds projected into the unc bases, re-relaxed), so that the repair is
   certified on the very atoms it must save — the formalization of the ad-hoc
   2026-11-07/08 Ti seeding that reproduced −848.40648875 exactly.
6. As the campaign, I want the 6-31G on-demand source, so that no repair depends
   on any cached campaign solution (the source converges fresh in seconds).
7. As the campaign, I want the projected seed threaded through `run_rohf` as a
   `dm0`, so that the seeded re-relax runs under the existing character-forced
   hooks with no new SCF path.
8. As the campaign, I want the density-direct quasi-mf variant certified, so that
   the `@s` source exists even if seeded SCFs still diverge.
9. As the campaign, I want `conditioned_basis` expressed as generalized
   contractions, so that the fallback cut is a standalone basis object the
   existing pipeline accepts unchanged.
10. As the physicist, I want the head floor measured (fresh heads on the cut
    basis vs the trusted anchors), so that the truncation's physical cost is a
    number, not a hope.
11. As the campaign, I want the twins' `@s` legs run and `_unc_s_block` retired,
    so that the tables complete and the guard's reason is gone.
12. As the campaign, I want the UHF anchor class re-tested under the seed, so
    that the same repair closes the second blocked row if it closes the first.
13. As the physicist, I want `subspace_rohf` (valence confined to the distilled
    span, core exact), so that each contraction's mean-field cost is a
    certificate rather than an inference from CI totals.
14. As the suite, I want the whole-AO-space degenerate anchor (the subspace ROHF
    equals the free ROHF exactly), so that the runner is pinned
    exact-by-construction at the one limit that admits it.
15. As the physicist, I want the per-contraction ΔE table (Sc first, Ti
    mirrored), so that the adoption decision un-gating the optimizer classes is
    measured.
16. As the physicist, I want the plan and its progress tracked in-repo, so that
    compaction loses nothing and any session resumes from the record.

## Implementation Decisions

**The seams and their homes** (deletion-test justified):

- `windows.project_orbitals(mol_b, mol_a, c_a, thresh=None) -> (c_b, report)` —
  the ONE primitive both repair uses lean on. Inside: the cross overlap
  `intor_cross('int1e_ovlp')`, the solve `S_B W = X C_A` **in the truncated
  eigenbasis** when `thresh` is given (the conditioning lives INSIDE the
  projector — `S_B⁻¹` amplifies the null space), Löwdin-orthonormalization of
  the projected span (`windows.lowdin`). `thresh=None` is deliberately the naive
  projector — the battery's stress-test foil. The report: per-column capture,
  the S_B spectrum at the cut, modes dropped — misuse visible, not silent.
  Home: `windows.py` beside `build_full_mo` (the S-metric column-set family;
  deleting it scatters the algebra into two real consumers — a real seam).
- `reference.projected_source(cfg, ref, source_basis="6-31g", thresh=<declared>)
  -> (c_occ_b, occ, report)` — the on-demand source adapter: forced-occupation
  ROHF on the small basis via the EXISTING `run_rohf` + `occupation_postcheck`
  (no basin gate there — the post-check is the honesty meter), then the
  projection. The threshold constant (Stage 1a's measured gap) is declared HERE
  — the one owner of the policy.
- The seed seam: `run_rohf(cfg, chkfile, ref=, dm0=)` → `_converged_cached(...,
  dm0=)` → `_accelerated` calling `m.run(dm0=dm0)` — three signatures, one
  thread, no new path. The re-relax stays under the character-forced hook
  (stronger than MOM — the hook IS the occupation anchor; no MOM machinery).
- `run_series(..., seed_basis=None)`; `CaseSpec.series_seed="6-31g"` declared on
  the twins — the repair is a declaration, not a code path.
- Use 2 — the density-direct quasi-mf (`mo_coeff`/`mo_occ` only) feeds
  `basis.valence_block` UNCHANGED (its only mf consumer). **FINALIZED**: stays
  sandbox-certified but unintegrated if the seeded branch wins on both twins
  (the virtmix pattern — architecturally possible, not declared).
- `basis.conditioned_basis(mol, thresh) -> (spec_dict, report)` — canonical
  orthogonalization per l-block as GENERALIZED CONTRACTIONS (one shell per l,
  all primitives as exponents, the retained eigenvectors as contraction rows) —
  a standalone raw PySCF spec, the `_SC_UNC_BASIS` kind, `build_mol`-ready
  today; the conditioned shells' labels are shell counters → the existing
  `"l0".."l9"` forced-occupation vocabulary. Home: `basis.py` (the
  contraction-recipe family).
  **AMENDED (Stage 1a, measured 2026-11-12 — the NO-GO branch): the twins'
  spectra show NO null space and NO gap** (every eigenvalue > 5e-4; l=2
  bottoms at 5e-2; the largest tail ratio jump ×3–4; `elements/<twin>/logs/
  spectrum.log`) — the cut design's premise fails at its own checkpoint, so
  **C is DEMOTED from a basis cut to the projector's internal guard**: the
  `thresh` solve stays (the robustness detail + the synthetic stress's foil)
  but removes nothing on these bases; a conservative default (1e-8-class) is
  the declared policy. **P is THE repair** (the seed first, the
  density-direct second); the divergence mechanism is the SCF dynamics on
  the amplified multi-basin landscape, not AO nullity. C's
  `conditioned_basis` is NOT built this wave (no input class exists for it
  now); the battery's synthetic near-singular anchors keep the projector's
  truncation certified.
- `reference.subspace_rohf(mf, c_valence, occ_closed, occ_open, ...)` — the
  `run_folded` mirror INVERTED: there the core is frozen and the valence solves
  free in the complement; here the VALENCE is confined to the distilled span
  while the CORE stays exact (its full AO family, free to relax — the cusp
  lesson). Mechanism: the same kernel-shim technique (ROHF unmodified on the
  S-orthonormal span `[_lane-core-family Löwdin'd ⊕ valence projected against
  it]`, `get_ovlp = I`, projected `h`, the veff shim), character-gated filling
  (`run_folded`'s pick-rule discipline). **FINALIZED DEFAULT**: lands as a
  SIBLING of `run_folded` sharing the tiny helpers (the bookkeepings genuinely
  differ — frozen-core field vs co-relaxing core); an extraction happens only if
  the deletion test shows the shared kernel is real at integration.

**The decision tree on the twins** (Stage 1b's output):

- P1 — the seeded forced SCFs converge on B → the real series on B; the `@s`
  legs run through `run_series` unchanged (best: full solutions in the complete
  span).
- P2 — P1 diverges → the density-direct quasi-mf (the `@s` source without any
  SCF on B).
- C — both fail → the conditioned basis; its solutions embed back into B.
- Certificates judge every draw (converged filling, `occupation_postcheck`, the
  term meters) — never the energy alone.

**Provenance** — no new tag grammar; the store payload records the source
family. **FINALIZED**: revisit only if the bench measures the constructible
variants' distilled radials disagreeing at the certificate floor (the
`s_r`-orthonormality/capture tolerances, `basis.py`).

**The sandbox protocol**: `sandbox/` in-repo, imported by nothing, deleted at
integration (git history keeps the bench). The filesystem rule honored — no
scratch outside the repo.

**Threading**: **64 OMP threads for the bench** (the user's directive
2026-11-12; the box carries 128 cores). The runs wave keeps OMP 16; a concurrent
bench + regen totals 80 threads — allowed, but the one-regen-at-a-time rule
stays.

## Testing Decisions

The repo's three-layer chain, mapped:

- (a) independent constructions — the projection cross-checked against the
  naive-vs-truncated pair on synthetic near-singular input;
- (b) external anchors — the trusted heads (both twins, more than one symmetry
  class);
- (c) exact analytic identities — the degenerate anchors below.

The battery:

| # | Anchor | Certifies |
|---|---|---|
| 1 | B = A | identity, exact — zero parameters |
| 2 | B ⊇ span(A) | the projected orbitals reproduce A's exactly |
| 3 | H → He → C | hand-checkable counts/energies |
| 4 | near-singular stress ± `thresh` | what the conditioning contributes; shared with C's synthetic near-duplicate bug-and-fix |
| 5 | the known-answer transfer | the contracted heads → the unc bases → re-relax → the trusted anchors |
| 6 | C's head floor | fresh heads on B′ vs trusted B anchors = the measured truncation cost |
| 7 | subspace: whole-AO-space == free ROHF; ΔE ≥ 0 monotone (r1 ⊂ s1 ⊂ s2); the closed-shell cross-check vs `confined_core_orbitals` | the runner exact at the degenerate limit, honest in between |

Tier placement: the synthetic anchors port to `tests/unit` (fast tier); the
twins' transfers and certificates to `tests/integration`. Every new construction
ships with a pin that would catch a sign/phase error. Per-commit gate: fast tier
+ the diff's targeted files (5–10 min); the full-suite run stays RETIRED; a red
is never the verdict — isolation reruns are.

## Out of Scope

- The runs wave itself (ticket 36 owns batches B→C→D→E→A; Stage 0 context only).
- The pinned-lane payload seam (map Issues; the declared-class decision pending).
- The same-l windows (`window_distilled`'s declared future seam — the
  member-preserving cross-shell orthogonalization).
- U14e's own session (the ECP record, separate per the handoff).
- The queued CASCI review (map Issues; called separately by the user).
- Spin–orbit, transition moments, g-factors (the candidate list beyond).

## Further Notes

Dependencies: **bench ⇒ the `@s` sources exist ⇒ the twins' `@s` legs run ⇒
`_unc_s_block` retires**; the subspace-ROHF exploration is independent of the
bench (different seams) but follows it in the declared order. The exploration's
outcome un-gates the two pending optimizer classes (the prior treatments'
re-homing, the twins' primary legs) — the anchor-first voids exist to wait for
it.

Costs: Stage 1a minutes; 1b ~an afternoon (battery + the twins' decision runs);
1c a 5–10 min gate; Stage 2 the D-batch class; Stage 3 minutes per solve.

Compaction: this ticket + `message.txt` are the record; the Progress section is
the live state — update it in the same commit as any state change.

## Progress

- [ ] **STAGE 0** (ticket 36; approved): batches B→C→D→E→A — re-render + commit
      per batch; U14e LAST.
- [x] **STAGE 1a** — the spectrum measurement: per-l-block S spectra on both
      twins, committed to `elements/<twin>/logs/spectrum.log`; **VERDICT =
      NO-GO on the clean-gap premise** (no null space: every eigenvalue
      > 5e-4, l=2 ≥ 5e-2; continuous decay, max jump ×3–4) — the cut
      refuted at its own checkpoint, C demoted to the projector's internal
      guard, P promoted to THE repair (see the AMENDED bullet above). The
      threshold policy: conservative (removes nothing here; the guard for
      genuinely singular inputs).
- [ ] **STAGE 1b** — the sandbox bench: `project_orbitals` + the battery
      (anchors 1–6) + the twins' decision runs (P1 → P2 → C).
- [ ] **STAGE 1c** — integration: the primitive/adapter/seed thread land
      (`windows` / `reference` / `basis` as won); the battery's synthetic
      anchors port to `tests/unit`, the transfers to `tests/integration`;
      lockstep docs; the sandbox deleted.
- [ ] **STAGE 2** — the `@s` wave: the twins' series declarations (seeded,
      censused e_refs) + `series_seed` + `_unc_s_block` retires; the CI anchors
      then the treatments; re-render + commit; the UHF anchor re-test.
- [ ] **STAGE 3** — the subspace-ROHF: `subspace_rohf` (sandbox → battery →
      integrate, anchor 7); the ΔE-vs-anchor certificate table (Sc first, Ti
      mirrored); the adoption decision recorded in map.md.
