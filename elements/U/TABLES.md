# U — uranium (f³d¹ outside the ECP core)

Stuttgart RLC ECP78 + cut-Dyall valence basis; 14 active electrons'
reference **ROHF −50.381852 Eh** (gated Newton; the reference
determinant is a numerically pure ⁵S° state — NOTES). Standard sector
M_S = 2. pur substrate mandatory (raw [H,L²] = 2.2; pur 5e-13).
Provenance: `results/u_14e/u_run*.log` (tiers, NO pilot, ticket-23
stages), `results/u_14e/u_sa*.log`, `u_six13.log`. NOTES.md carries
the ticket-23 narrative and open questions.

## T3 — Window × treatment (energies in Eh; ΔE vs frozen anchors in
mEh; ‖·‖ = many-body commutator; |[P,·]| = span certificate; labels
from exact eigh in the converged window — PySCF ci vectors straggle
in the 17-fold ⁵L cluster)

| window | dets | \|[P,Lz]\| | \|[P,L⊥]\| | CASCI, Eh | ⟨L²⟩ | relax | SA (17×⁵L) |
|---|---|---|---|---|---|---|---|
| 1f+1d band | 495 | ~1e-13 | ~1e-13 | −49.459 | 42.000 | — | — |
| 2f+1d | 3,876 | ~1e-13 | ~1e-13 | −49.751 | 72.000 | — | — |
| 2f+2d | 10,626 | ~1e-13 | ~1e-13 | −49.894 | 72.000 | — | — |
| 3f+2d | 31,465 | ~1e-13 | ~1e-13 | −50.081 | 72.000 | — | — |
| 4f+2d | 73,815 | ~1e-13 | ~1e-13 | **−50.392** | 72.000 | — | — |
| ROHF MO (open+virtuals, sizes 12–100) | 495+ | O(7–9) at every size | O(7–9) | −50.382…−50.384 | 0.000–0.258 | — | — |
| NO groups 12 col (d5+f7) | 495 | ~1e-13 | ~1e-13 | −49.889731 | 72.000 | −49.899729 (−10.0 mEh), ‖·‖ 1.3e-13 | **−50.423810, ⟨L²⟩ 72.0000 exact** (fold ‖·‖ 1.3e-4) |
| NO groups 17 col | 2,380 | ~1e-13 | ~1e-13 | −49.890855 | 72.000 | −49.900366 (−9.5 mEh) | degraded: span opens 4.3/5.4, fold ‖·‖ 10 |
| NO groups 24 col (= full window) | 10,626 | ~1e-13 | ~1e-13 | −49.894287 | 72.000 | — | — |
| NO-12 + 7s (CAS(6,13), (5,1)) | 16,731 | 2e-13 | 4e-13 | −49.894866 | 72.000 | — | — |

Free single-root CASSCF: broken at every start and size tried (NO
5/12/17/24, shell 12/24) — window lost in one macro step, gradient
stuck ~1e-2, never converged, ⟨L²⟩ 0.6–31. See NOTES.

Descent race: Δf pays f³d¹ −339/−330 mEh (at 1d/2d) but f²d² only
−187; Δd pays f²d² −144 but f³d¹ only −69/−60. At 4f+2d the ⁵L° (f³d¹)
manifold stands below the best f²d² root by > 37 mEh — ground term
**⁵L° (f³d¹, L = 8, ⟨L²⟩ = 72)**.
