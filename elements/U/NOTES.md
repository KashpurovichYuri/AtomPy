# U — NOTES (factual, dated; details behind TABLES.md)

- **The ⁵S° reference (open; tickets 21/24).** The gated ROHF
  determinant is a numerically pure ⁵S° state: ⟨L²⟩ = 8e-10, ⟨L_z⟩ = 0,
  parity −1, in three independent constructions; all 3,876 generic
  dets of the 1f+1d window sit at ⟨L²⟩ ≥ 6. "Below ROHF" = beating an
  ⁵S° mean field. Mechanism open (average-of-configuration smearing
  hypothesis). Landscape support: free optimization descends THROUGH
  the ⁵S° region (macro step 1 of the NO-12 free run lands at
  −50.3828, the ROHF energy).
- **Ticket 23 — relaxation vs growth (2026-10-24).** SA over the
  17-root ⁵L block from the NO-12 start: **−50.423810 at exactly
  quantized ⁵L** (exact-eigh label; fold ‖[H,Lz]‖ 1.3e-4; span cert
  8e-4) — 42 mEh below the ⁵S° ROHF reference; the first U object
  below the mean field at certified symmetry. Mechanism measured as
  RADIAL: window l-content exactly preserved (l2 = 5.0, l3 = 7.0,
  l0 = l1 = 0.000) while 3.53 columns of weight move to f/d primitives
  outside the selected exponent bands (band weight 7.000 → 4.706) —
  ticket-24's missing-radial axis recovered by relaxation, at 495 dets.
  Reproduced bit-exactly across independent reruns.
- **Free rungs die on U.** From every start: window identity lost in
  macro step 1 (−493 mEh), orbital gradient stuck at ~1e-2, energy
  crawling ~−18 µEh/cycle at 100 cycles (never converged), ⟨L²⟩ broken
  (0.6–31), span certs O(5–8), ‖[H,Lz]‖ up to 270. Same-span starts
  (shell vs NO) land 0.14 mEh apart with chaotic gauges. The ⁵S°
  attractor takes the money and destroys the terms; only the ensemble
  holds the gauge (SA-12). SA-17 overreaches — pulls 5 more f members
  in (l3: 7 → 12, still l-pure) but m-incomplete: certs 4.3/5.4.
- **Relax rungs.** NO-12: −10.0 mEh vs frozen at machine-exact closure;
  NO-17: −9.5. Core relaxation ≈ 2× the entire frozen 12→24 growth
  axis (+4.6 mEh). Truncated runs chaotic across reruns; converged
  endpoints bit-stable.
- **7s² is inert as valence (six13).** CAS(6,13) exact eigh:
  −49.894866 vs campaign −49.894287 (−0.58 mEh), ⟨N7s⟩ = 1.97,
  ⟨Nf⟩ = 2.00, ⟨Nd⟩ = 2.03, span machine-exact; ground unchanged
  f²d² ⁵L. The accidental 6p-member promotion (identification bug,
  fixed) moved −551 mEh with an m-incomplete lone member: INNER core
  pairs are the next correlation tier, to be done m-complete.
- **NO groups (ticket 22).** Term-average of the 17 ⁵L roots → groups
  d/f/d/f = 5/7/5/7, occupations 0.399/0.286/5e-4/2e-4; group windows
  exactly closed at every size; 12 cols within 4.6 mEh of the full
  window where the AO lattice at equal dets sits +431 mEh higher.
- **Label rule (hard-won).** PySCF ci vectors straggle in the 17-fold
  degenerate cluster: mc.ci mislabeled ⟨L²⟩ = 27.4 where the exact
  ground in the same window/fold is 72.0000. Diagonal-population
  expectations are |c_i|²-weighted, never signed sums.
- **Open.** ⁵S° mechanism; inner core-pair valence (m-complete
  design); f³d¹-vs-f²d² ordering inside one NO window; NIST; 3rd-d
  point; the window × treatment product (ticket 20).
