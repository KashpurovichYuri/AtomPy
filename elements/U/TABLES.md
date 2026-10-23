# U — the uranium campaign (u_14e)

Stuttgart RLC ECP78 + cut-Dyall valence basis; 14 electrons outside the
ECP core; ground reference **ROHF −50.381852 Eh** (gated Newton; the
reference determinant is a numerically pure ⁵S° state, ⟨L²⟩ = 0 in
three independent constructions — mechanism open, AOC-smearing
hypothesis). Standard sector M_S = 2 (the cheapest quintet). pur-core
substrate mandatory (raw [H, L²] = 2.2; pur 5.1e-13).

Provenance: `results/u_14e/u_run_part2.log` (tiers), `u_run_no.log`
(NO pilot), `u_sa_casscf.log` / `u_sa_diag.log` / `u_sa12.log` /
`u_sa17.log` / `u_sa_relax.log` (ticket 23), `u_six13.log` (7s² probe).
Cross-atom verdicts: README §7.5.

## Window lattice — exact symmetry vs raw MOs at matched cost

| window | CAS | dets | AO, Eh | ⟨L²⟩ | ROHF MO, Eh | ⟨L²⟩ | NO |
|---|---|---|---|---|---|---|---|
| 1f+1d ‖ open+8 virt. | (4,12) | 495 | −49.459 | 42.000 | −50.382 | 0.000 | −49.890 |
| 2f+1d ‖ open+15 | (4,19) | 3,876 | −49.751 | 72.000 | −50.383 | 0.132 | — |
| 1f+3d ‖ open+18 | (4,22) | 7,315 | −49.733 | 42.000 | −50.384 | 0.192 | — |
| 2f+2d ‖ open+20 | (4,24) | 10,626 | −49.894 | 72.000 | −50.384 | 0.254 | −49.894 |
| 3f+1d ‖ open+22 | (4,26) | 14,950 | −49.993 | 72.000 | −50.384 | 0.258 | — |
| 3f+2d | (4,31) | 31,465 | −50.081 | 72.000 | — | — | — |
| 4f+1d | (4,33) | 40,920 | −50.332 | 72.000 | — | — | — |
| 4f+2d | (4,38) | 73,815 | **−50.392** | 72.000 | — | — | — |

Descent race: Δf pays f³d¹ −339/−330 mEh (at 1d/2d) but f²d² only
−187; Δd pays f²d² −144 but f³d¹ only −69/−60. At 4f+2d the ⁵L° (f³d¹)
manifold stands below the best f²d² root by > 37 mEh — the ground term
is **⁵L° (f³d¹, L = 8, ⟨L²⟩ = 72)**, odd parity.

Raw-ROHF MO spans never rotation-close (‖[P,L⊥]‖ = 7–9 at every size
4..100); the reference determinant is representable only from 4f+2d up
(3f+2d projection norm 0.948).

## Natural-orbital groups (ticket 22)

Term-average of the 17 ⁵L roots (2f+2d window) → four occupation groups
d/f/d/f = 5/7/5/7, occupations 0.399 / 0.286 / 0.0005 / 0.0002. Group
windows are exactly rotation-closed at every size (certs ~1e-13):
12 columns (495 dets) hold ⁵L at −49.890 (+4.6 mEh vs full window); the
AO lattice at the same 495 dets sits +431 mEh higher.

## Optimizer rungs — relaxation vs growth (ticket 23)

All labels from exact eigh in the converged windows (PySCF ci vectors
straggle in the 17-fold ⁵L cluster: mc.ci mislabels ⟨L²⟩ as 27.4 where
the exact ground is 72.0000).

| rung | window | CAS | dets | E, Eh | ⟨L²⟩ | symmetry |
|---|---|---|---|---|---|---|
| frozen CASCI (anchor) | NO d5+f7 | (4,12) | 495 | −49.889731 | 72.0000 | machine-exact |
| frozen growth | all 4 NO groups | (4,24) | 10,626 | −49.894287 | 72.0000 | machine-exact |
| core relax (span frozen) | NO-12, core relaxed | (4,12) | 495 | −49.899729 | 72.0000 | machine-exact |
| SA-CASSCF, 17×⁵L | NO-12 start, radials relaxed | (4,12) | 495 | −50.423810 | 72.0000 | fold 1.3e-4, span 8e-4 |
| 7s² promoted | NO-12 + 7s pur column | (6,13) | 16,731 | −49.894866 | 72.0000 | machine-exact |
| free single-root | any start (NO nk=1..3, shell 12/24) | — | — | ≈ −50.43 | broken | unconverged; spans O(5–8) |

Verdicts:

- **The deep correlation is radial.** SA keeps the window exactly in
  the f/d shells (l-content l2 = 5.0, l3 = 7.0 unchanged, l0 = l1 =
  0.000) while 3.53 columns of weight move to f/d primitives OUTSIDE
  the selected exponent bands — ticket-24's missing-radial axis
  recovered by relaxation instead of growth: −534 mEh at exactly
  quantized ⁵L. First U object below the ⁵S° ROHF mean field at
  certified symmetry (−42 mEh).
- **Only the ensemble survives the ⁵S° attractor.** Free single-root
  loses the window in one macro step (−493 mEh), never converges
  (gradient ~1e-2 stuck), lands broken (⟨L²⟩ 0.6–31). SA-17
  overreaches: pulls 5 more f members in (l3: 7 → 12, still l-pure)
  but the span OPENS (4.3/5.4) and the fold breaks ([H,Lz] = 10.2,
  ⟨L²⟩ = 71.88) — l-pure is not m-complete.
- **Core relaxation** at frozen span: +10.0 mEh at machine-exact
  closure — ~2× the entire 12→24 frozen growth axis (+4.6 mEh).
  Truncated runs are chaotic across reruns; converged endpoints
  reproduce bit-exactly.
- **7s²-as-valence is inert**: −0.58 mEh, ⟨N7s⟩ = 1.97, ⟨Nf⟩ = 2.00,
  ⟨Nd⟩ = 2.03; the ground stays the f²d² ⁵L block, span machine-exact.
  (An accidental 6p-member promotion instead moved −551 mEh with an
  m-incomplete lone member — indicative that INNER core pairs are the
  next correlation tier, to be done m-complete.)

## Open questions

- The ⁵S° reference mechanism (why the gated ROHF det is an exact ⁵S°
  state; AOC-smearing hypothesis) — the free-run attractor is
  landscape-level support that ⁵S° is real on U.
- Inner core-pair valence (6p/…), done m-complete: the −551 mEh class.
- f³d¹-vs-f²d² ordering inside ONE NO window (ticket-22 remainder).
- NIST comparison; the 3rd-d U point.
