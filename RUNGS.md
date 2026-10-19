# The orbital-treatment ladder, in plain words

*(The derivation-level narrative is the README technical note — §5 locates
each rung's mechanism, §7 the measured numbers.)*

The "ladder" is one series of calculations on the same atom. The active
window (the few orbitals we correlate — e.g. carbon's three 2p) starts
identical in every rung; what differs is **what we allow to happen to the
orbitals**, from "touch nothing" to "optimize everything". Comparing the
rungs then tells us which orbital freedom produced which effect — energy
gain, symmetry loss, or nothing at all. Each row of the ladder report is
one rung.

The symmetry meter throughout: `[H, Lz]` and `[H, L²]` commutator norms.
They are exactly zero when the one-particle space (core + window) is a
complete set of angular-momentum functions; every rung's deviation is a
measured number, not a judgement.

## The CASCI rungs (window fixed, only the core field differs)

**casci_raw — "the honest SCF core".** We diagonalize in the raw shell
AOs and fold in the ROHF core field *as converged*. That core is
secretly lopsided: the SCF saw one specific open-shell electron (say px),
so the core density carries a tiny quadrupole-shaped polarization.
Result: commutators ≠ 0 at a small level and the term multiplet splits
in the `{M_L=0}` vs `{M_L=±1}` pattern. *Question answered: how big is
the raw reference's symmetry breaking?* (C: 1e-5 Eh; Sc: percent-level.)

**casci_avg — "the spherical version of the same core".** Identical, but
the core density is spherically averaged (its l=0 component) before
folding. Commutators drop to ~1e-13. *Question: was the breaking entirely
the core's anisotropy?* (Yes, on C/Si/Sc.)

**casci_frac — "a different, intrinsically spherical reference".** A
fractional-occupation SCF (each open-shell orbital occupied by 2/3 of an
electron, say) is spherical *by construction* — no averaging needed. But
it is a different physical reference (valence charge smeared vs
localized), so its energy differs at the mEh scale. *Question: how much
does the choice of reference itself matter?*

## The CASSCF rungs (the window may change; how freely?)

**casscf_relax — "window frozen, core may relax".** Only core↔virtual
rotations are allowed: the core polarizes in response to the correlated
active density, no orbital ever enters or leaves the window. *Question:
how much does plain core relaxation contribute, and at what symmetry
cost?*

**casscf_radial — "radial mixing only".** A window slot may rotate only
with same-(l, m) partners: 2px may borrow 3px, never py, never an s.
Every slot keeps its angular identity; the active electrons can only
spread radially (in–out correlation). *Question: how much correlation is
radially reachable with zero angular freedom?* (On C/Si: all of it.)

**casscf_virtmix — "no core exchange, virtuals free".** The window may
borrow any virtual character but never trade with the core. *Question:
is the Sc escape reachable through virtuals alone?* (Yes — blocking core
exchange does not prevent the escape; blocking l-change does.)

**casscf_shell / casscf_pyscf — "free optimization"** from the raw-shell
start / PySCF's default (canonical) start. *Question: where does an
unconstrained optimizer actually go, and does the starting guess matter?*
(On Sc it leaves the intended window entirely — the "escape".)

## The SA rungs (average over the term's M_L components)

**sa — "state-averaged, free".** The optimizer minimizes the average
energy of the 2L+1 multiplet components, which makes the density a
rotational scalar and the symmetric orbitals a stationary point.
*Question: does averaging restore the symmetry?*

**sa_radial — "state-averaged + radial mask".** Identity and member
structure frozen by the mask; SA only picks the (symmetric) gauge.
*Question: with everything else frozen, does SA drive both commutators to
zero?* (C: yes, 8e-11. Sc: down to the input-window baseline, ~1e-3.)

**sa_pyscf — "state-averaged from the canonical start".** *Question: is
SA's outcome start-dependent?* (Very much so on Sc.)

## The purified-substrate rungs (optional: `run_ladder(..., purified=True)`)

The core every CASSCF rung must orthogonalize its window against is NOT
rotationally closed (fact 17): the canonical core MOs carry small
m-asymmetric tails, so the *input window* already carries a closure
defect (Si 3.4e-3, Sc 2.1e-2) that every optimized rung inherits as a
commutator baseline. The `*_pur` rungs redo the same treatments on a
**purified core** — the canonical core replaced by the closest exactly
rotation-closed span (one shared radial per closed shell, replicated on
every m member; `purified_core_orbitals`). The physical reference is the
same scalar core the analysis fold already uses; the purification costs
only the discarded core polarization (C +3 µEh, Si +0.13 mEh,
Sc +0.67 mEh).

**casci_avg_pur / casscf_relax_pur / sa_relax_pur — "exact by
construction or by frozen window".** *Question: with the baseline gone,
which optimized rungs close at machine precision?* (Answer: all of the
frozen-window ones, on every atom — including **Sc SA-CASSCF with
unbroken Lz and L² at 1e-13** via sa_relax_pur.)

**sa_pur / sa_radial_pur — "free SA / radial SA on the closed
substrate".** (C/Si: machine-exact — the radial direction has a
restoring force, so SA's stationary point IS the symmetric window.
Sc: the member gauge is exactly flat — m-symmetric radial rotation gains
nothing, fact 18b — so the optimizer stops at a ~1e-6 tilt that no
convergence tightening removes; an optimizer floor, not a baseline. The
free-sa 4s↔3p swap attractor survives the purification: it is an
optimizer basin, not a core-span artifact.)

## Reading the report

- `E_ground`, `dE_vs_avg` — ground energy and its change vs casci_avg.
- `<L^2>` — angular momentum of the ground state (2 → triplet P for C).
- `|[H,Lz]|`, `|[H,L2]|` — the symmetry meters (see above).
- `win_l` — the window's l-purity (1.000 = no s/d leakage; necessary but
  NOT sufficient for symmetry — the m-members must also stay complete).
- `|[P,Lz]|`, `|[P,Lp]|` — the one-particle closure defects of the window
  span itself (see the README glossary appendix): what the many-body
  commutators can and
  cannot see, with no core fold and no CI. Both ~0 = the window span is
  certified rotation-closed.
- `aniso` — the core density's discarded anisotropy (casci rungs).
- `max|dE|` — our reconstructed energies vs PySCF's own: a bookkeeping
  check with known physical offsets (a polarized core we average but
  PySCF uses raw shows up here).
- `ortho` — orthonormality of the rung's MO matrix. Above 1e-6 the rung
  is flagged **DIVERGED** and its numbers are meaningless (this is the
  ticket-09 guard; it has a dedicated negative-control test).

## Known reading hazards

- Energies of CASSCF-family rungs vs casci rungs carry the
  window-definition offset (fact 14): compare within families.
- Single-run commutators of *free* single-root rungs are chaotic — they
  measure which member gauge the optimizer happened to land in (fact 16).
- An L≈0 state has small `<L^2>` in ANY window — only the commutators
  measure closure (the Sc escape "masquerades" as symmetric).
