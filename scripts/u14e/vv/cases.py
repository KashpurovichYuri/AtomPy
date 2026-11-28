"""The U14e CAS ladder the table carries.

A :class:`Case` is one table tag: how its active space is built and
what the committed records certify about it.  The tags and their
groups are ``vv.table``'s; this module gives each tag its window.
The core width and the sector come from the construction -- the
key's for a key window, the core's for a visible one -- never from a
hand-typed constant.

Every window is built on the trusted reference's purified core (the
core the optimization starts from; ``core_leak`` in each solution's
certificate reads how far the optimizer moved it).  The plain
``(4,12)`` and ``(4,15)`` tags take the VISIBLE windows (the shells'
own AO columns); every other tag takes a KEY window (the production
recipe on a cell address), except the ``@r1`` (6,13) composition,
whose s channel is replaced by the reference's own closed pair.

THE MINIMAL FAMILY (``vv.table.LADDER_MIN``, the window-variant
sweep) is built otherwise: its base is the free state average's
CONVERGED ``(4,12)@r1`` active block, assembled on the key window's
purified core -- the ground term's basin.  The frozen key window
cannot be that base: its lowest 17 roots are the 5L term of
5f^2 6d^2 (2.0 electrons in f, 2.0 in d; the ``@r1`` radial is the
mean over the configuration series), and a state average seeded there
tracks that term.  A tag ``(4,N)m+<channels>`` names the channels set
beside the base (``channels_of``): ``s`` an s virtual (8s; 9s lands
the same window), ``p`` a p radial, ``d`` a second d radial, ``f`` a
second f radial -- the window of each size with the lowest
state-averaged energy in the sweep.  The family's start carries the
free optimizer's residual asymmetry (rotation-closed to ~7e-04); the
constrained route's lm-pure basis projects it out.

Layering
--------
``TERM`` -- the declared term; :class:`Case`; ``CASES`` (the ladder)
and ``select_cases``; the constructors ``visible_window``,
``window_6p_r1`` and the minimal family's ``channels_of`` /
``minimal_window`` / ``converged_base`` (``window_from_key`` and
``window_6_13_r1`` are ``vv.windows``').
"""

import re
from dataclasses import dataclass
from functools import partial
from typing import Callable

import numpy as np

from atomterm.basis import series_window
from atomterm.core import purified_core_orbitals
from atomterm.windows import build_full_mo, occupied_core_idx

from vv.channels import CHANNELS, channel_block
from vv.reference import Substrate
from vv.routes import free_sa
from vv.table import LADDER_MIN
from vv.windows import (ActiveSpace, Built, composition_lines,
                        window_6_13_r1, window_from_key, window_from_shells,
                        window_with_channels)

TERM = (8, "u")               # the declared term: L = 8, ungerade
NROOTS = 2 * TERM[0] + 1      # its M_L components: the state average's roots

# eigenvalues the expectations name: <L^2> = L(L+1), <S^2> = S(S+1)
L2_TERM, L2_L9 = float(TERM[0] * (TERM[0] + 1)), 90.0
S2_QUINTET, S2_SEPTET = 6.0, 12.0


@dataclass(frozen=True)
class Case:
    """One table tag: its window and the records' expectations.

    ``window`` -- ``(substrate) -> Built``.

    ``frozen`` -- the frozen ground's certified (<L^2>, <S^2>) on the
        free route (None: any definite L).
    ``sa_term`` -- (<L^2>, <S^2>) every free-SA root must carry
        (None: the roots partition a flat manifold draw-dependently;
        disclosed).
    ``anchor`` -- the committed payload cell the free SA's root-0 E(2)
        is checked against (the declared-roots discipline).
    ``pinned_frozen`` / ``pinned_sa`` -- <L^2> every pinned component
        must carry (the pin-on-pur gate); None: no gate.
    ``pinned_class`` -- the pinned SA's component classes as the lane
        record certified them: (count within 0.01 of L2_L9, count
        within 0.01 of L2_TERM).
    ``s_pair`` -- the frozen ground must doubly occupy the s channel
        (the (6,13)@r1 acceptance control).
    ``frozen_below`` -- the frozen ground's dE against the reference
        must stay below this many mEh (a window that freezes far
        above its reference has abandoned the intended term).
    ``ground_config`` -- (electrons in the base's f block, in its d
        block) the ensemble must hold: the configuration is the
        ground one (5f^3 6d^1), which (<L^2>, <S^2>) cannot tell from
        another of the same term.
    ``radial`` -- the constrained state average runs on the window
        (``radial_term``: the expectations every one of its roots
        carries); ``pinned`` -- the pinned routes run on it.
    ``start_closed`` -- the start window is a constructed,
        rotation-closed one (the closure gate); False where the start
        carries an optimizer's residual.
    ``casscf_cycles`` -- the single-root free CASSCF's macro-cycle
        budget (None: the library's); a case whose single-root run
        goes to the cap declares a shorter one.
    """

    tag: str
    window: Callable[[Substrate], Built]
    frozen: tuple[float, float] | None = None
    sa_term: tuple[float, float] | None = None
    anchor: str | None = None
    pinned_frozen: float | None = None
    pinned_sa: float | None = None
    pinned_class: tuple[int, int] | None = None
    s_pair: bool = False
    frozen_below: float | None = None
    ground_config: tuple[float, float] | None = None
    radial: bool = False
    radial_term: tuple[float, float] | None = None
    pinned: bool = True
    start_closed: bool = True
    casscf_cycles: int | None = None


# ------------------------------------------------------------ the windows

SHELLS_412 = (("10f", 7), ("13d", 5))
SHELLS_415 = SHELLS_412 + (("14p", 3),)
SECTOR_4 = (4, 0)             # the 5f3 6d1 configuration, Ms = 2

ADDR = {
    "(4,16)": "u14e_10f13d14p7s_4-0_casci_ao_-_pur_avg",
    "(6,13)": "u14e_10f13d7s_5-1_sa-casscf_ao_-_pur_avg",
    "(4,12)@r1": "u14e_5f6d@r1_4-0_casci_ao_-_pur_avg",
    "(4,15)@r1": "u14e_5f6d14p@r1_4-0_casci_ao_-_pur_avg",
    "(4,16)@r1": "u14e_5f6d14p7s@r1_4-0_casci_ao_-_pur_avg",
    "(6,13)@r1": "u14e_5f6d7s@r1_5-1_casci_ao_-_pur_avg",
}


def core_pairs(sub: Substrate, sector) -> int:
    """The core pairs that leave ``sector`` active -- derived from
    the molecule's electron count, never hand-typed (a hand-typed
    count once shifted every meter slice by one column)."""
    n = sub.mol.nelectron - sum(sector)
    if n <= 0 or n % 2:
        raise ValueError(f"sector {sector} leaves {n} electrons for "
                         f"the core -- not a positive even count")
    return n // 2


def visible_window(sub: Substrate, shells, sector) -> Built:
    """A visible window on the purified core: the family's own
    construction (:func:`vv.windows.window_from_shells`)."""
    ncore = core_pairs(sub, sector)
    core = purified_core_orbitals(sub.mol, sub.mf, ncore)[0]
    return window_from_shells(sub.mol, sub.mf, shells, core, sector)


# ------------------------------------------------- the minimal family

F_BLOCK, D_BLOCK = 7, 5     # the base's columns: one 5f radial, one 6d radial
# a tag's channel tokens -> the declared channels they name
MIN_TOKENS = {"s": ("8s",), "2s": ("8s", "9s"), "p": ("7p",),
              "d": ("7d",), "f": ("5f",)}
_BASE_TAG = "(4,12)@r1"
_BASE: dict = {}            # the converged base, per substrate (in-process)


def channels_of(tag: str) -> tuple[str, ...]:
    """The channels a minimal-family tag names beside the base.

    Raises
    ------
    ValueError : not a minimal-family tag; an unknown token; the
        tag's size disagrees with its channels' widths.
    """
    m = re.fullmatch(r"\(4,(\d+)\)m(?:\+(.*))?", tag)
    if not m:
        raise ValueError(f"{tag!r} is not a minimal-family tag "
                         f"(the form is (4,N)m+<channels>)")
    tokens = m.group(2).split("+") if m.group(2) else []
    unknown = [t for t in tokens if t not in MIN_TOKENS]
    if unknown:
        raise ValueError(f"unknown channel token(s) {unknown} in {tag!r}; "
                         f"declared: {sorted(MIN_TOKENS)}")
    channels = tuple(c for t in tokens for c in MIN_TOKENS[t])
    width = F_BLOCK + D_BLOCK + sum(CHANNELS[c].l * 2 + 1 for c in channels)
    if width != int(m.group(1)):
        raise ValueError(f"{tag!r}: its channels make {width} active "
                         f"orbitals, not the tag's {m.group(1)}")
    return channels


def converged_base(sub: Substrate) -> Built:
    """The minimal family's base: the free state average's CONVERGED
    ``(4,12)@r1`` active block, assembled on the key window's
    purified core (the module docstring) -- computed once per
    substrate.

    Raises
    ------
    ValueError : the free state average's certificate refuses.
    """
    key = id(sub)
    if key not in _BASE:
        frozen = CASES[_BASE_TAG].window(sub)
        sol = free_sa(sub, frozen.space, NROOTS)
        st = sol.states[0]
        mo = build_full_mo(sub.mf, occupied_core_idx(sub.mf, st.ncore),
                           np.array(st.active), c_core=frozen.space.core)
        _BASE[key] = Built(
            ActiveSpace(mo, st.ncore, st.ncas, st.nelecas),
            (f"base: the free state average's converged {_BASE_TAG} active "
             f"block on its purified core (ensemble {sol.ensemble:.8f}, "
             f"{sol.conv})",))
    return _BASE[key]


def minimal_window(sub: Substrate, channels) -> Built:
    """A minimal-family window: the channels beside the converged
    base (:func:`vv.windows.window_with_channels`)."""
    base = converged_base(sub)
    return window_with_channels(
        sub, base, [channel_block(sub, CHANNELS[c]) for c in channels])


def window_6p_r1(sub: Substrate, k: int) -> Built:
    """The 5f/6d/6p composition with per-shell k from the ``@r1``
    density blocks -- a composition the key grammar cannot carry (its
    k is uniform per window), built by a direct distillation call.

    Raises
    ------
    ValueError : the tail certificate refuses a k beyond the radial
        signal (the declared refusal class).
    """
    mol, mf, spec = sub.mol, sub.mf, sub.spec
    ncore = core_pairs(sub, SECTOR_4)
    bands = {sh: spec.shell_bands[sh] for sh in ("5f", "6d", "6p")}
    act, rep = series_window(mol, [mf], mf.get_ovlp(), ["5f", "6d", "6p"],
                             bands, {"5f": 1, "6d": 1, "6p": k}, ncore)
    core = purified_core_orbitals(mol, mf, ncore)[0]
    mo = build_full_mo(mf, occupied_core_idx(mf, ncore), act, c_core=core)
    return Built(ActiveSpace(mo, ncore, act.shape[1], SECTOR_4),
                 tuple(composition_lines(rep)))


# -------------------------------------------------------------- the ladder

CASES: dict[str, Case] = {c.tag: c for c in (
    Case("(4,12)", partial(visible_window, shells=SHELLS_412, sector=SECTOR_4),
         frozen=(42.0, S2_QUINTET), sa_term=(L2_TERM, S2_QUINTET),
         anchor="u14e_10f13d_4-0_sa-casscf_ao_-_pur_avg"),
    Case("(4,15)", partial(visible_window, shells=SHELLS_415, sector=SECTOR_4),
         frozen=(42.0, S2_QUINTET),
         anchor="u14e_10f13d14p_4-0_sa-casscf_ao_-_pur_avg"),
    Case("(4,16)", partial(window_from_key, address=ADDR["(4,16)"]),
         frozen=(42.0, S2_QUINTET)),
    Case("(6,13)", partial(window_from_key, address=ADDR["(6,13)"]),
         frozen=(12.0, S2_SEPTET), pinned_class=(16, 1)),
    Case("(4,12)@r1", partial(window_from_key, address=ADDR["(4,12)@r1"]),
         frozen=(L2_TERM, S2_QUINTET), sa_term=(L2_TERM, S2_QUINTET),
         pinned_frozen=L2_TERM, pinned_sa=L2_TERM),
    Case("(4,15)@r1", partial(window_from_key, address=ADDR["(4,15)@r1"]),
         frozen=(L2_TERM, S2_QUINTET)),
    Case("(4,16)@r1", partial(window_from_key, address=ADDR["(4,16)@r1"]),
         frozen=(L2_TERM, S2_QUINTET)),
    Case("(4,21)@r1", partial(window_6p_r1, k=3)),
    Case("(6,13)@r1", partial(window_6_13_r1, address=ADDR["(6,13)@r1"]),
         frozen=(L2_TERM, S2_QUINTET), s_pair=True, frozen_below=400.0),
    # the minimal family: the ground term is the frozen ground, the
    # constrained state average keeps it, the ensemble is 5f^3 6d^1;
    # its single-root free CASSCF runs to the cap (observed 200 cycles
    # at <L^2> 64.7 -- the root leaves the term), on a 50-cycle budget
    *(Case(tag, partial(minimal_window, channels=channels_of(tag)),
           frozen=(L2_TERM, S2_QUINTET), ground_config=(3.0, 1.0),
           radial=True, radial_term=(L2_TERM, S2_QUINTET), pinned=False,
           start_closed=False, casscf_cycles=50)
      for tag in LADDER_MIN),
)}


def select_cases(tags=()) -> list[Case]:
    """The declared cases in table order, or the named ones (an
    unknown tag is refused)."""
    from vv.table import table_keys
    unknown = [t for t in tags if t not in CASES]
    if unknown:
        raise ValueError(f"unknown case tag(s) {unknown}; the ladder is "
                         f"{list(CASES)}")
    order = [k.split("|")[-1] for k in table_keys()]
    wanted = order if not tags else tags
    return [CASES[t] for t in dict.fromkeys(wanted) if t in CASES]


__all__ = ["ADDR", "CASES", "Case", "D_BLOCK", "F_BLOCK", "L2_L9",
           "L2_TERM", "MIN_TOKENS", "NROOTS", "SECTOR_4", "S2_QUINTET",
           "S2_SEPTET", "SHELLS_412", "SHELLS_415", "TERM", "channels_of",
           "converged_base", "core_pairs", "minimal_window",
           "select_cases", "visible_window", "window_6p_r1"]
