"""The channel map (the correlation-channel exploration)."""

import numpy as np  # noqa: E402

# the channel map (the correlation-channel exploration) ---------------

def band_chains_by_suffix(mol):
    """The AO index groups per primitive band and per harmonic type:
    returns (bands, types) -- ``bands[band_label] = [ao indices]``
    (the label's numeric+letter prefix is the primitive band; "10f-2"
    and "10f+2" share the band "10f"), ``types[band_label] =
    {type_suffix: [ao indices]}`` (the harmonic type within the band:
    the m-type identification -- "pz"/"py"/"px", "dz^2"/"dxy"/...,
    the f's signed m: "f+0"/"f-2"/...)."""
    import re as _re
    labels = [lab.split()[-1] for lab in mol.ao_labels()]
    bands = {}
    types = {}
    for i, lab in enumerate(labels):
        m = _re.match(r"(\d+[spdf])(.*)$", lab)
        band, suffix = m.group(1), m.group(2)
        bands.setdefault(band, []).append(i)
        types.setdefault(band, {}).setdefault(suffix, []).append(i)
    return bands, types


def band_population(mol, mf, bands, weights=None):
    """Per MO, the population fraction on each band (the S-metric
    band weights of the column, normalized over the bands).
    ``weights`` -- the MO occupation weights (default: the
    occupations).  Returns (band_order, (n_mo, n_bands) fractions)."""
    s = mf.get_ovlp()
    c = np.asarray(mf.mo_coeff)
    order = sorted(bands)
    mat = np.zeros((c.shape[1], len(order)))
    for k, band in enumerate(order):
        idx = bands[band]
        sub = s[np.ix_(idx, idx)]
        cs = c[idx, :]
        mat[:, k] = np.einsum("ai,ab,bi->i", cs, sub, cs).real
    tot = mat.sum(axis=1, keepdims=True)
    tot[tot == 0] = 1.0
    if weights is not None:
        w = np.asarray(weights, dtype=float)
        if w.ndim == 1:
            mat = mat * w[:, None]
    return order, mat / tot


def family_spectrum(mol, mf, shell, all_bands, ncol, families_types):
    """The occupied-signal Gram spectrum over ONE shell family's full
    primitive set (the device's own m-summed valence density --
    ``valence_block`` through the family's type chains), plus the
    per-band attribution of each radial eigenvector: the radial's OWN
    coordinates expanded to the slots (the raw |c|^2 per band,
    normalized).  The attribution is deliberately NOT the S-metric:
    the strongly-overlapping same-l primitives spread the S-metric
    weight uniformly (~1/n each -- uninformative), while the raw
    coordinates resolve the radial's primitive composition.  Returns
    (eigenvalues descending, the per-band weights per eigenvector:
    (n_eig, n_bands))."""
    from atomterm.basis import valence_block
    # the chains: ONE group per harmonic type, across ALL the family
    # bands' slots of that type (the valence_block sums the type
    # blocks -- each (n_primitives x n_primitives) -- giving the
    # family's radial spectrum)
    type_slots = {}
    type_bands = {}
    for band in all_bands:
        for suffix, idx in sorted(families_types[band].items()):
            type_slots.setdefault(suffix, []).extend(idx)
            type_bands.setdefault(suffix, []).extend([band] * len(idx))
    type_order = sorted(type_slots)
    chain = [type_slots[t] for t in type_order]
    slot_list = [(i, band, t) for t, ty in enumerate(type_order)
                 for i, band in zip(type_slots[ty], type_bands[ty])]
    block = valence_block(mf, chain, ncol)
    lam, tvecs = np.linalg.eigh(block)
    lam = lam[::-1].real
    tvecs = tvecs[:, ::-1].real
    ubands = sorted(set(b for _, b, _ in slot_list))
    fam_idx = [i for i, _, _ in slot_list]
    s = mf.get_ovlp()
    s_fam = s[np.ix_(fam_idx, fam_idx)]
    wbands = np.zeros((len(lam), len(ubands)))
    for e, tvec in enumerate(tvecs.T):
        # the radial, expanded from the type basis to the family slots
        v_fam = np.zeros(len(slot_list))
        for n, (_, _, t) in enumerate(slot_list):
            v_fam[n] = tvecs[t, e]
        # the per-band attribution: the eigenvector's OWN coordinates
        # (the raw |c|^2) -- the S-metric would spread uniformly over
        # the strongly-overlapping same-l primitives (measured ~1/n)
        for k, band in enumerate(ubands):
            bidx = [n for n, (_, b, _) in enumerate(slot_list)
                    if b == band]
            wbands[e, k] = float(np.sum(v_fam[bidx] ** 2))
    tot = wbands.sum(axis=1, keepdims=True)
    tot[tot == 0] = 1.0
    wbands /= tot
    return lam, (ubands, wbands)
