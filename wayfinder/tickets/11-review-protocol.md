# Ticket 11: review protocol (/code-review, /improve-codebase-architecture)

Status: open (run after the atoms wave 06–09 lands)

## Question

How do we review the accumulated wave and find deepening opportunities?

## Resolution (protocol)

Skills already installed in `.opencode/skill/`: `code-review`
(model-invoked, registered), `improve-codebase-architecture`
(user-invoked; follow its SKILL.md + HTML-REPORT.md manually).

**Step 1 — /code-review** (after tickets 06–09 merged):

- Fixed point: commit `507b796` (review the atoms wave only; the earlier
  wave 0a36d56..507b796 was continuously test-validated during pairing).
- Two axes, run as PARALLEL sub-agents (per the skill):
  - *Standards*: AGENTS.md conventions (module table accuracy, hard-won
    facts, docstring/comment discipline, test naming) + Fowler smell
    baseline (duplication, long parameter lists, shallow modules).
  - *Spec*: faithfulness to wayfinder tickets 06–09 (every promised
    behavior present, no scope creep beyond tickets).
- Findings → fix immediately or file as tickets; re-run until clean.

**Step 2 — /improve-codebase-architecture** (after review fixes):

- Survey for deepening opportunities. Anticipated candidates (verify, don't
  assume): shared parts of detbasis/secondquant (the parts/tensor layer),
  drivers.py scope growth (orbital selection + purification + SCF + CAS),
  pipeline residues (run_carbon vs run_atom), test duplication between
  test_operators.py and test_atoms.py.
- Produce the visual HTML report per the skill; present candidates to the
  user; GRILL through the chosen ones before changing anything (the skill
  is a survey, not a rescue).
- Agreed changes become new wayfinder tickets (one session each).

**Step 3**: record the protocol (fixed point convention, cadence: review
after each ticket wave) in AGENTS.md so future sessions apply it.
