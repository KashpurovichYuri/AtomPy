# Ticket 42 — `perturb.py` rungs 0–1: the FOIS and Epstein–Nesbet

Spec: `wayfinder/specs/pt2-perturbation.md` (the wave's contract —
the vocabulary: EN, MP, Dyall, FOIS, intruder, cumulants,
size-consistency — lives in the module docstring at implementation,
general form, no atom names).

## Goal

The module born: the coupling-only FOIS enumeration (rung 0) and the
Epstein–Nesbet E(2) in two independent constructions — the
sum-over-states formula and the first-order linear solve (rung 1,
name-symmetric twins). H₀ = diag(H) in the external space; the frozen
core excluded (core orbitals never excited). The interface: one
function per rung + the shared `PT2Result` (correction, meters:
floor = min|E⁰ − H₀μμ|, w_ref = 1/(1 + Σ|C⁽¹⁾μ|²), FOIS completeness
vs the brute-force count).

## Status: READY (nothing blocks)

## Inputs

- `detbasis.py` — the Slater–Condon machinery; the coupling column's
  twin is a brute-force `slater_condon_matrix` slice at test sizes;
- `operators.py` — `full_space_integrals`, the substrate signature,
  the ci_stack ordering bridge;
- `tests/carbon_reference.py` + the `carbon` fixture — the frozen
  substrate (NO SCF in tests);
- synthetic substrates for the degenerate anchors (random
  non-symmetric integrals — every term class exercised).

## Verification

(a) twins elementwise: EN sum vs EN solve; FOIS coupling column vs
the brute-force H-matrix slice. (c) degenerate anchors: CAS = whole
space → FOIS empty → E(2) = 0 exactly; one active electron → the
closed answer; the empty core. Cardinality vs the brute-force
excitation count at pinned sizes. The carbon fixture pins the first
real EN number (the observed fresh value in a comment). The fast
tier + `tests/unit/test_perturb.py` = the per-commit gate.
