# Ticket 15: review protocol (/code-review, /improve-codebase-architecture)

Status: open (runs FIRST in the reordered wave: before 14 and Ti)

## Question

How do we review the accumulated wave and find deepening opportunities?

## Resolution (protocol)

Skills already installed in `.opencode/skill/`: `code-review`
(model-invoked, registered), `improve-codebase-architecture`
(user-invoked; follow its SKILL.md + HTML-REPORT.md manually).

**Step 1 — /code-review** (review everything through ticket 11):

- Fixed point: commit `f5cfb6f` (review the tickets 10-11 delta; the
  earlier waves were continuously test-validated during pairing).
- Two axes, run as PARALLEL sub-agents (per the skill):
  - *Standards*: AGENTS.md conventions (module table accuracy, hard-won
    facts, docstring/comment discipline, test naming) + Fowler smell
    baseline (duplication, long parameter lists, shallow modules).
  - *Spec*: faithfulness to wayfinder tickets 10-11 (every promised
    behavior present, no scope creep beyond tickets).
- Findings -> fix immediately or file as tickets; re-run until clean.
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
