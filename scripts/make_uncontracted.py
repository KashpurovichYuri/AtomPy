"""Split a contracted basis into its single-primitive blocks.

The uncontracted variant of a basis file: every contraction block
becomes one block per primitive, coefficient 1.0, the l-section and
primitive order preserved. The output is the nwchem format PySCF's
``gto.parse`` reads (the ``data/u_dyall_cut_basis.txt`` precedent).

The split is MECHANICAL and loses the contraction grouping by design
-- the chemical shells of the split basis are declared separately (the
shell layer, ``CaseSpec.shell_bands``) and certified by the measured
composition (``scripts/band_composition.py``): on an uncontracted
basis the AO labels are radial counters, not shells.

Usage:
    python scripts/make_uncontracted.py Sc ccpvdz > data/sc_ccpvdz_unc.txt

The script verifies its own output before printing: re-parsing the
emitted text must reproduce the primitive exponents exactly, block
for block, every block single-coefficient.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pyscf import gto  # noqa: E402

_ANG = "spdfgh"


def uncontracted_blocks(basis: list):
    """The basis's blocks as single-primitive blocks, order preserved.

    Parameters
    ----------
    basis : a PySCF basis list (``gto.basis.load`` output):
        general-contraction blocks ``[l, [exp, coef1..coefK], ...]``
        (K contraction coefficients per primitive).

    Returns
    -------
    ``[[l, [(exp, 1.0)], ...]`` -- one block per primitive.
    """
    out = []
    for block in basis:
        l = block[0]
        # the general-contraction rows: [exp, coef1, ..., coefK] -- the
        # K contraction coefficients are DISCARDED by the split (one
        # block per primitive, coefficient 1.0)
        for row in block[1:]:
            out.append([l, [(float(row[0]), 1.0)]])
    return out


def render_nwchem(symbol: str, blocks: list) -> str:
    """The blocks as nwchem-format text (one ``<El> <l>`` header per
    block, exponent and coefficient indented)."""
    lines = []
    for l, prims in blocks:
        lines.append(f"{symbol} {_ANG[l].upper()}")
        for exp, coef in prims:
            lines.append(f"    {exp:.10E}  {coef:.1f}")
        lines.append("")
    return "\n".join(lines)


def main(argv) -> int:
    if len(argv) != 3:
        print(__doc__.strip().splitlines()[-3])
        return 2
    symbol, source = argv[1], argv[2]
    basis = gto.basis.load(source, symbol)
    blocks = uncontracted_blocks(basis)
    text = render_nwchem(symbol, blocks)

    # the self-check: the emitted text re-parses to the same
    # primitives, every block single-coefficient
    reparsed = gto.parse(text)
    assert [b[0] for b in reparsed] == [b[0] for b in blocks]
    assert all(len(b[1:]) == 1 for b in reparsed)
    assert [b[1][0] for b in reparsed] == [b[1][0][0] for b in blocks]

    sys.stdout.write(text)
    print(f"# uncontracted split of {symbol}/{source}: "
          f"{len(blocks)} primitive blocks", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
