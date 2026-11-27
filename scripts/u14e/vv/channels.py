"""HF-virtual correlation channels: a shell of the reference's VIRTUAL
orbitals, promoted into an active space beside a base window as one
more radial of its (l, member) set.

A :class:`Channel` names the reference's virtual columns that carry
a shell and the AO bands realizing it.  :func:`channel_block` turns
it into ``2l + 1`` S-orthonormal columns; :func:`window_with_channels`
(``vv.windows``) sets them beside a key window.  The construction is
measured, never assumed: :func:`virtual_spectrum` reads each virtual
orbital's angular character, purity and extent.

TWO CONSTRUCTIONS, by angular momentum:
  l = 0 -- the virtual orbital's s SECTOR: a single-member shell is
      its own radial.  (The raw orbital of a C1 reference carries a
      small d admixture -- observed weights d 3.4e-04 and 4.7e-04 --
      and the closure meter is first order in that amplitude: 6e-02
      for the raw orbital, 2e-14 for its s sector.  The distillation
      device's radial is unreliable on the s family of an uncontracted
      basis: its row-block Gram overshoots where the diffuse functions
      are nearly dependent -- observed captures 0.28 and 0.04 of the
      two lowest s virtuals.)
  l >= 1 -- the distillation device (``windows.window_distilled``) on
      the shell's l-pure members: ONE m-averaged radial, replicated on
      the full m-set -- rotation-closed by construction whatever the
      anisotropy or the splitting of the members (a split shell's
      missing members simply do not feed the signal).

Layering
--------
:class:`Channel` and ``CHANNELS`` (the declared inventory) ;
:func:`virtual_spectrum` ; :func:`channel_block`.
"""

from dataclasses import dataclass

import numpy as np

from atomterm.angular import l2_matrix
from atomterm.windows import (l_per_row, shell_weights, window_distilled)

from vv.reference import Substrate

_L = "spdfgh"


@dataclass(frozen=True)
class Channel:
    """One channel of the reference's virtual space.

    ``shell`` -- its label; the last letter is the angular momentum.
    ``columns`` -- the reference's virtual MO columns carrying it.
    ``bands`` -- the AO bands realizing it (the device's layer).
    ``tail`` -- the radial certificate: the discarded Gram spectrum
        must sit below this fraction of the leading eigenvalue (a
        declared relaxation carries its measured value).
    ``pure`` -- the minimum l-weight a column must carry to feed the
        radial signal (a strongly mixed column is left out and
        reported; a weakly mixed one feeds the l rows only).
    """

    shell: str
    columns: tuple[int, ...]
    bands: tuple[str, ...]
    tail: float = 1e-4
    pure: float = 0.97

    @property
    def l(self) -> int:
        return _L.index(self.shell[-1])


_P = ("16p", "17p", "18p", "19p")
_S = ("16s", "17s", "18s", "19s", "20s")
_D = ("12d", "13d", "14d", "15d")
_F = ("9f", "10f", "11f", "12f", "13f", "14f")

# The declared inventory of the U14e reference (the ROHF rung's
# virtuals by angular character; energies relative to the HOMO, a.u.
# of <r^2>).  Each tail cites the measured second-to-leading Gram
# ratio of its shell.
CHANNELS: dict[str, Channel] = {c.shell: c for c in (
    Channel("7p", (0, 1, 2), _P, tail=0.01),              # +0.028 Eh x3; <r^2> 103..111; ratio 1.2e-03
    Channel("8s", (3,), _S),                             # +0.052 Eh; <r^2> 114
    Channel("7d", (4, 5, 9, 10, 11), _D, tail=0.2),       # +0.070..0.122 Eh x5 (z^2 apart); ratio 1.4e-01
    Channel("8p", (6, 7, 8), _P, tail=0.01),              # +0.073..0.077 Eh x3; ratio 7.3e-04
    Channel("8d", (12, 13, 14, 15), _D, tail=0.01),       # +0.147..0.154 Eh x4; ratio 1.5e-03
    Channel("9s", (16,), _S),                            # +0.191 Eh; <r^2> 75
    Channel("6f", (17, 18, 19, 20, 23, 25, 26, 27), _F,   # +0.206..0.219 Eh; ratio 4.2e-02
            tail=0.05),                                     # (6 members feed; 2 are p-mixed)
    Channel("5f", (28, 29, 30, 31), _F, tail=0.01),       # +0.284..0.291 Eh x4 (split); <r^2> 12..13; ratio 3.8e-04
)}


def virtual_spectrum(sub: Substrate, below: float = 0.9) -> list[dict]:
    """The reference's virtual orbitals by angular character.

    Parameters
    ----------
    sub : the reference substrate.
    below : list the virtuals up to this many hartree above the HOMO.

    Returns
    -------
    One record per virtual, in energy order: ``index`` (the MO
    column), ``de`` (hartree above the HOMO), ``l`` (the dominant
    angular momentum's letter), ``weight`` (its S-metric l-weight),
    ``l2`` (<L^2>; l(l+1) for a pure orbital), ``r2`` (<r^2>).
    """
    mol, mf = sub.mol, sub.mf
    s = mf.get_ovlp()
    c = np.asarray(mf.mo_coeff)
    occ = np.asarray(mf.mo_occ)
    energy = np.asarray(mf.mo_energy)
    homo = float(np.max(energy[occ > 0]))
    idx = np.flatnonzero((occ == 0) & (energy - homo < below))
    idx = idx[np.argsort(energy[idx])]
    names, w = shell_weights(mol, c[:, idx], s)
    l2, r2 = l2_matrix(mol), mol.intor("int1e_r2")
    out = []
    for k, j in enumerate(idx):
        col = c[:, j]
        norm = float(col @ s @ col)
        li = int(np.argmax(w[k]))
        out.append({"index": int(j), "de": float(energy[j] - homo),
                    "l": names[li], "weight": float(w[k, li]),
                    "l2": float(np.real(col @ l2 @ col)) / norm,   # L^2 is stored complex-typed
                    "r2": float(col @ r2 @ col) / norm})
    return out


def channel_block(sub: Substrate, channel: Channel) -> tuple[np.ndarray, list[str]]:
    """The channel's ``2l + 1`` S-orthonormal columns (see the module
    docstring for the two constructions).

    Returns
    -------
    (columns, report) -- (nao, 2l+1) and the construction's
    disclosure lines.

    Raises
    ------
    ValueError : no column is l-pure enough to carry the shell; or
        the device's radial certificate refuses (propagated).
    """
    mol, mf = sub.mol, sub.mf
    s = mf.get_ovlp()
    c = np.asarray(mf.mo_coeff)[:, list(channel.columns)]
    names, w = shell_weights(mol, c, s)
    letter = _L[channel.l]
    keep = [k for k in range(c.shape[1]) if w[k, channel.l] >= channel.pure]
    mixed = [channel.columns[k] for k in range(c.shape[1]) if k not in keep]
    if not keep:
        raise ValueError(
            f"channel {channel.shell!r}: none of its columns "
            f"{list(channel.columns)} carries >={channel.pure:g} "
            f"{letter}-weight (best {w[:, channel.l].max():.3f})")
    lines = [f"channel {channel.shell}: {len(keep)} l-pure signal "
             f"column(s) {sorted(channel.columns[k] for k in keep)}"
             + (f"; left out as {letter}-mixed: {mixed}" if mixed else "")]
    if channel.l == 0:
        # the s sector of each orbital (an l = 0 direction is
        # L-invariant; the C1 gauge's d tail is projected out)
        block = c[:, keep] * (l_per_row(mol) == 0)[:, None]
        kept = np.einsum("ai,ab,bi->i", block, s, block).real
        block = block / np.sqrt(kept)
        lines.append(
            "channel: the virtuals' s sector (a single-member shell is "
            "its own radial; the C1 gauge's d tail is projected out -- "
            "keeps cos^2 " + " ".join(f"{x:.5f}" for x in kept) + ")")
        return block, lines
    block, rep = window_distilled(mol, c[:, keep], s, [channel.shell],
                                  {channel.shell: channel.bands},
                                  {channel.shell: 1}, tail=channel.tail)
    r = rep[channel.shell]
    lam = r["spectrum"]
    lines.append(
        f"channel: one m-averaged radial on {len(channel.bands)} bands; "
        f"the Gram's ratio 2nd/1st {lam[1] / lam[0]:.1e} (declared tail "
        f"{channel.tail:g}); capture "
        + " ".join(f"{x:.2f}" for x in r["capture"]))
    return block, lines


__all__ = ["CHANNELS", "Channel", "channel_block", "virtual_spectrum"]
