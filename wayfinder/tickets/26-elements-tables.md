# Ticket 26: elements/ tables to the approved T3a/T3b formats

Status: OPEN (after ticket 25 -- the cross-check column's meaning
depends on it)

Origin: the 2026-10-24 review. The user reinstated the AO/MO/NO lattice
comparison as the headline table ("the difference and importance of
proper window choice") and rejected the earlier merged window x
treatment table that buried it.

## T3a -- window lattice at matched cost (frozen CASCI)

Header line: the ROHF reference (value + configuration + gated note).
Columns: window | CAS | dets | AO, Eh | <L^2> | MO, Eh | <L^2> |
NO, Eh | <L^2> — <L^2> is a SEPARATE column after every E column.
Caption: AO = shells/primitive bands; MO = the ROHF open MOs +
virtuals ranked by valence-band weight; NO = occupation groups of the
term-averaged density. <L^2> exact on AO/NO at every size, broken on
MO — the window choice IS the symmetry statement.

Approved Ti example (values verified against the campaign logs):

ROHF reference −848.406489 Eh (3d²4s², gated)

| window | CAS | dets | AO, Eh | ⟨L²⟩ | MO, Eh | ⟨L²⟩ | NO, Eh | ⟨L²⟩ |
|---|---|---|---|---|---|---|---|---|
| 3d | (4,5) | 50 | −848.054 | 6.000 | — | — | — | — |
| 3d4s | (4,6) | 120 | −848.395 | 12.000 | −848.405 | 0.000 | −848.399 | 12.000 |
| 3d4d | (4,10) | 1,200 | −848.232 | 6.000 | — | — | — | — |
| 3d4s4p | (4,9) | 756 | −848.429 | 12.000 | — | — | −848.432 | 12.000 |
| 3d4s4d | (4,11) | 1,815 | −848.409 | 12.000 | −848.414 | 0.010 | — | — |
| 3d4s4p4d | (4,14) | 5,096 | −848.438 | 12.000 | −848.420 | 0.201 | −848.438 | 12.000 |
| 3d4s4p4d5d | (4,19) | 18,411 | −848.439 | 12.000 | — | — | — | — |
| 3d4s4p4d4f | (4,21) | 27,930 | −848.446 | 12.000 | — | — | — | — |

## T3b -- treatment ladder at the standard window

ONE energy convention (as run = PySCF's own variational total; the
comparable column; ΔE = E − ROHF). "Our fold" energies are NOT a table
column (analysis ruler only; its offset is the cross-check column).
Columns: treatment | E (as run), Eh | ΔE (mEh) | ⟨L²⟩ | ‖[H,Lz]‖ |
‖[H,L²]‖ | ‖[P,Lz]‖ | ‖[P,L⊥]‖ | cross-check (mEh).

Approved Ti example (values verified against
results/Ti_ccpvdz_s2_ladder.json):

| treatment | E (as run), Eh | ΔE (mEh) | ⟨L²⟩ | ‖[H,Lz]‖ | ‖[H,L²]‖ | ‖[P,Lz]‖ | ‖[P,L⊥]‖ | cross-check (mEh) |
|---|---|---|---|---|---|---|---|---|
| ROHF (anchor) | −848.406489 | 0 | — | — | — | — | — | 0 |
| CASCI, averaged | −848.395940 | +10.5 | 12.0000 | 2.2e-14 | 4.1e-13 | 5.9e-15 | 7.6e-15 | 89.5 |
| CASCI, raw | −848.395940 | +10.5 | 12.0000 | 4.9e-02 | 2.6e-01 | 5.9e-15 | 7.6e-15 | 88.2 |
| CASCI, fractional | −848.358519 | +48.0 | 12.0000 | 1.2e-13 | 7.5e-13 | 5.9e-15 | 7.6e-15 | 84.9 |
| CASSCF, free | −848.437740 | −31.3 | 0.0503 | 4.3e-01 | 1.1e+00 | 2.8 | 4.0 | 2.8 |
| CASSCF, core-relax | −848.396109 | +10.4 | 12.0 | 2.3e-08 | 1.7e-07 | 3.2e-08 | 5.9e-08 | 2.7 |
| SA-CASSCF, free | −848.416563 | −10.1 | 12.003 | 1.1e+00 | 6.0e+00 | 1.4 | 1.5 | 851 |
| SA-CASSCF, core-relax, pur | −848.394753 | +11.7 | 12.0 | 5.2e-13 | 2.9e-12 | 7.2e-15 | 1.2e-14 | 0.0 |

Legend: E (as run) = the run's own variational total; ΔE vs the ROHF
anchor; cross-check = |E(our analysis fold) − E(as run)| max over
roots — near zero only on the pur substrate; tens of mEh where the raw
core polarizes; ~850 mEh where the partition moved (the swap — totals
only, §5-L4). Rows whose window identity changed (escape, swap) are
different problems — in as-run totals the variational ordering is
sane. ⟨L²⟩ on the raw fold staying 12.0000 while the commutators jump
is §2.4 in one number (expectations shift at ε²).

## Scope

All five atoms (C, Si, Sc, Ti, U): TABLES.md rewritten to T1 / T2 /
T3a / T3b (+ U: the descent-race note and the ticket-23 table in T3b
style; C/Si: T3b as rung lists as before but in the new column
format). Sweep markers (<!-- study:begin/end -->) added for the future
product rows. README §7.3 adopts the T3b format for the mechanism
tables.

## Values provenance

Every number from the campaign artifacts (results/*_ladder.json +
logs); the Ti T3b verified in-session 2026-10-24; others re-verified
at execution.
