# Ticket 15: review protocol (/code-review, /improve-codebase-architecture)

Status: open (runs FIRST in the reordered wave: before 14 and Ti)

## Question

How do we review the accumulated wave and find deepening opportunities?

## Resolution (protocol)

Skills already installed in `.opencode/skill/`: `code-review`
(model-invoked, registered), `improve-codebase-architecture`
(user-invoked; follow its SKILL.md + HTML-REPORT.md manually).

**Step 1 — /code-review** (review everything through ticket 11): DONE

- Fixed point `f5cfb6f`; both axes run as parallel sub-agents.
- Standards: 0 hard violations; fixed in the follow-up commit: invariant
  naming (`test_quantum_number_round_trips_measured_momenta`,
  `test_build_mol_yields_origin_centered_spherical_atom`), dead import,
  `_s_ortho_err` signature, /tmp chk paths -> session tmp fixtures,
  SC_EXACT_TERMS shared via conftest, same-basin helper, the `if name ==
  "C"` fixture special-case dropped (run_ladder for all cases).
- Spec: the loosened-assertion finding above (fixed, disclosed in ticket
  10); the missing S^2-spectrum and [Ar]-core pins added
  (`s_values` case key; `test_sc_rohf_reference_is_3d1_4s2_with_argon_core`);
  README/.gitignore/wave-reorder confirmed user-directed, not creep.
- New physics-coverage test: `test_sc_frac_occ_selects_open_shell_by_projection`
  pins BOTH the projection-selected {3d,4s} block and the energy-ordered
  failure mode (p-containing ensemble) -- the fact-15(c) mechanism.

- Fixed point: commit `f5cfb6f` (review the tickets 10-11 delta; the
  earlier waves were continuously test-validated during pairing).
- Two axes, run as PARALLEL sub-agents (per the skill):
  - *Standards*: AGENTS.md conventions (module table accuracy, hard-won
    facts, docstring/comment discipline, test naming) + Fowler smell
    baseline (duplication, long parameter lists, shallow modules).
  - *Spec*: faithfulness to wayfinder tickets 10-11 (every promised
    behavior present, no scope creep beyond tickets).
- Findings -> fix immediately or file as tickets; re-run until clean.
  (Known findings folded in pre-emptively: the "2S-like" wording — now
  "dominantly-L=0 basin of the escaped window, optimizer artifact" — and
  the README refresh; the relax-rung comm attribution remains OPEN, filed
  under ticket 14 pre-items.)
  Known findings to fold in: the "2S-like" wording (say "dominantly-L=0
  basin of the escaped window, optimizer artifact, not a physical term"),
  the relax-rung comm (~1e-3) attribution, README staleness.

**Step 2 — /improve-codebase-architecture** (after review fixes):

- Survey for deepening opportunities. Anticipated candidates (verify, don't
  assume): shared parts of detbasis/secondquant (the parts/tensor layer),
  drivers.py scope growth (orbital selection + frac-occ + SCF + CAS +
  masks), the pipeline/ladder fold-and-analyze duplication, run_ladder
  orchestration accretion, test duplication between test_operators.py and
  test_atoms.py.
- Produce the visual HTML report per the skill; present candidates to the
  user; GRILL through the chosen ones before changing anything (the skill
  is a survey, not a rescue).
- Agreed changes become new wayfinder tickets (one session each). The
  ticket-14 pre-items (rung registry, constraint spec, dedupe) are prime
  candidates -- the survey should confirm or reshape them.

**Step 3**: record the protocol (fixed point convention, cadence: review
after each ticket wave) in AGENTS.md so future sessions apply it.
