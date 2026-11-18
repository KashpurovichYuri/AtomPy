# Ticket 45 — `PERTURBATION_TABLE.md`: the declared-scope render

Spec: `wayfinder/specs/pt2-perturbation.md` (the schema section owns
the table grammar).

## Goal

`elements/{C,Si}/perturbation.json` + `PERTURBATION_TABLE.md`: ONE
correction table, raw columns only — `window | sector | core |
root (term) | ndet | E_ref | ⟨L²⟩ | variant | E_PT2 | E_total |
ΔE_total | floor | w_ref | notes` — rendered under the DECLARED
scope (default: SA-CASSCF roots with identifiable terms across the
window ladder; the m-incomplete MO/MOC controls excluded; everything
else joinable by declaration). The README reading key (R, the
guards, the conditional verdicts) lands in the same commit — the
mechanisms' one home.

## Status: BLOCKED BY 43 (44's variant columns join at its close;
the schema declares all variant columns from day one)

## Inputs

- `perturbation.json` (the 42/43/44 outputs, keyed window × sector ×
  core × root × variant);
- `scripts/build_methods_tables.py` — the emitter conventions to
  mirror (byte-determinism pins, the marker pattern).

## Verification

Renderer pins: byte-determinism, declared-scope policy (the
lane/lattice emitter pattern). The tables render from validated
rungs only (ticket 43's green gate). The fast tier stays green.
