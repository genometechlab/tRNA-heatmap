"""
Heatmap generation and pileup data export for tRNA pileup data.

Public API
----------
  plot(arrays, sprinzl_axis, ref_to_sprinzl, output_path, **kwargs)
      Generate heatmap from raw pileup arrays + Sprinzl mapping.

  plot_from_npz(npz_path, output_path, **kwargs)
      Generate heatmap directly from a pre-computed NPZ file.

  delta_heatmap(npz_paths, output_prefix, **kwargs)
      Compute pairwise delta heatmaps across a list of NPZ files.

  delta_from_arrays(arrays_by_stem, sprinzl_axis, ref_to_sprinzl, output_prefix, **kwargs)
      Compute pairwise delta heatmaps directly from pileup arrays (no NPZ required).

  plot_sprinzl_coverage(sprinzl_axis, ref_to_sprinzl, output_path, **kwargs)
      Presence/absence heatmap showing which Sprinzl positions each tRNA covers.
      No pileup data required — driven by cmalign output alone.

  save_pileup(arrays, sprinzl_axis, ref_to_sprinzl, base_path)
      Save pileup data as both <base_path>.npz and <base_path>.tsv.

Shared keyword arguments (plot / plot_from_npz / delta_heatmap / delta_from_arrays)
---------------------------------------------------------------
  palette          : str   'light-high' (default) | 'dark-high' | 'viridis'
  ylabel           : str   Y-axis label. Default 'tRNA Reference Names'.
  title            : str   Figure title. Default 'tRNA Alignment Pileup Heatmap'.
  show_insertions  : bool  Show insertion columns (Xi1, Xi2 ...) in the plot.
                           They are always saved in NPZ/TSV. Default False.

Cell rendering
--------------
  Colored cell   : position has coverage; color = mismatch rate.
  Blank cell     : position exists in the tRNA but had zero reads (no coverage).
  Black dot (●)  : cmalign alignment shows no base at this Sprinzl position for
                   this tRNA. In delta plots, a dot appears if EITHER tRNA in the
                   pair lacks the position.
"""

import itertools
import os
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['pdf.fonttype'] = 42   # embed fonts as TrueType in PDF (editable in Illustrator)
plt.switch_backend('agg')           # non-interactive backend; safe on headless servers


# ---------------------------------------------------------------------------
# Palette registry
# ---------------------------------------------------------------------------

PALETTES = {
    'light-high': 'YlGnBu',   # yellow=low, dark blue=high mismatch (default)
    'dark-high':  'hot_r',     # white=0 %, black=100 % mismatch
    'viridis':    'viridis',   # perceptually uniform, purple → yellow
}

_DELTA_CMAP = 'Spectral_r'    # always used for delta plots (diverging)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _plot_aligned(matrix, ref_names, sprinzl_axis, output_path,
                  palette='light-high',
                  ylabel='tRNA Reference Names',
                  title='tRNA Alignment Pileup Heatmap',
                  show_insertions=False,
                  vmin=0, vmax=1,
                  cmap_override=None,
                  no_base_mask=None,
                  cbar_label=None,
                  dpi=300):
    """
    Render and save a heatmap from a pre-aligned metric matrix.

    Parameters
    ----------
    matrix : np.ndarray, shape (n_refs, n_sprinzl)
        Metric values (e.g. mismatch rate or delta). NaN = no coverage or no base.
    ref_names : list[str]
    sprinzl_axis : list[str]
    output_path : str
    palette : str
        One of PALETTES keys. Ignored when cmap_override is set.
    ylabel : str
    title : str
    show_insertions : bool
        If False (default), columns whose Sprinzl label contains 'i<digit>'
        (e.g. '36i1') are stripped before rendering. They remain in NPZ/TSV.
    vmin, vmax : float
        Color scale limits. Use vmin=-1, vmax=1 for delta plots.
    cmap_override : str or None
        If set, overrides palette (used for delta plots).
    no_base_mask : np.ndarray (bool), shape (n_refs, n_sprinzl), or None
        True where the tRNA has no base at that Sprinzl position (cmalign gap).
        Black dots are drawn only at these cells. If None (old NPZ files), falls
        back to drawing dots at all NaN cells.
    cbar_label : str or None
        Override the colorbar axis label. Defaults to 'Mismatch rate' or
        '\u0394 Mismatch rate' based on vmin.
    """
    # --- Filter insertion columns ---
    if not show_insertions:
        keep = [i for i, lbl in enumerate(sprinzl_axis)
                if not re.search(r'i\d+', lbl)]
        sprinzl_axis = [sprinzl_axis[i] for i in keep]
        matrix = matrix[:, keep]
        if no_base_mask is not None:
            no_base_mask = no_base_mask[:, keep]

    n_refs, n_cols = matrix.shape
    cmap_name = cmap_override if cmap_override else PALETTES.get(palette, 'hot')
    cmap_obj = plt.get_cmap(cmap_name).copy()
    cmap_obj.set_bad(alpha=0)   # NaN cells: fully transparent

    # --- Figure sizing: aim for roughly square cells ---
    cell_size = 0.25                          # inches per cell
    fig_w = max(6.0, n_cols * cell_size + 3.0)
    fig_h = max(3.0, n_refs * cell_size + 2.0)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # pcolormesh renders vector rectangles in PDF (no raster interpolation blur)
    im = ax.pcolormesh(matrix, cmap=cmap_obj, vmin=vmin, vmax=vmax)
    ax.set_xlim(0, n_cols)
    ax.set_ylim(n_refs, 0)   # top-to-bottom row order, matching imshow convention

    # --- Black dots: no-base positions only (or all NaN for old files) ---
    if no_base_mask is not None:
        dot_rows, dot_cols = np.where(no_base_mask)
    else:
        dot_rows, dot_cols = np.where(np.isnan(matrix))

    if dot_rows.size > 0:
        # pcolormesh cell centres are at col+0.5, row+0.5
        ax.scatter(dot_cols + 0.5, dot_rows + 0.5, s=20, c='black', marker='o',
                   zorder=2, linewidths=0)

    cbar = fig.colorbar(im, ax=ax, shrink=0.5)
    cbar.ax.tick_params(labelsize=6)
    if cbar_label is None:
        cbar_label = 'Mismatch rate' if vmin == 0 else '\u0394 Mismatch rate'
    cbar.ax.set_ylabel(cbar_label, fontsize=6, labelpad=10, rotation=270, va='bottom')

    # pcolormesh cell centres at i+0.5
    ax.set_yticks([i + 0.5 for i in range(n_refs)])
    ax.set_yticklabels(ref_names, fontsize=6)
    ax.tick_params(axis='y', length=4)

    ax.set_xticks([i + 0.5 for i in range(n_cols)])
    ax.set_xticklabels(sprinzl_axis, rotation=90, fontsize=7)
    ax.tick_params(axis='x', length=6)

    ax.set_xlabel('Sprinzl position')
    ax.set_ylabel(ylabel)
    ax.set_title(title, pad=7)

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, transparent=True)
    plt.close(fig)
    print(f"Heatmap saved to {output_path}")


def _build_aligned(arrays, sprinzl_axis, ref_to_sprinzl):
    """
    Map pileup arrays to Sprinzl coordinate columns.

    Returns
    -------
    aligned : np.ndarray, shape (4, n_refs, n_sprinzl)
        axis 0 = [match, mismatch, insertion, deletion]
    ref_names : list[str]
        Sorted reference names (row order in aligned).
    no_base_mask : np.ndarray (bool), shape (n_refs, n_sprinzl)
        True where the tRNA has no base at that Sprinzl position.
    """
    ref_names     = sorted(arrays.keys())
    n_refs        = len(ref_names)
    n_sprinzl     = len(sprinzl_axis)
    sprinzl_index = {label: i for i, label in enumerate(sprinzl_axis)}
    aligned       = np.zeros((4, n_refs, n_sprinzl), dtype=np.float64)
    no_base_mask  = np.ones((n_refs, n_sprinzl), dtype=bool)

    for row_i, name in enumerate(ref_names):
        arr = arrays[name]
        for ref_pos, label in enumerate(ref_to_sprinzl.get(name, [])):
            col = sprinzl_index.get(label)
            if col is not None and ref_pos < arr.shape[1]:
                aligned[:, row_i, col]  = arr[:, ref_pos]
                no_base_mask[row_i, col] = False

    return aligned, ref_names, no_base_mask


def _mismatch_rate(entry, common_refs, common_axis):
    """
    Extract a mismatch-rate matrix aligned to common_refs × common_axis.

    entry : dict with keys 'data' (4,n,m), 'refs' list[str], 'axis' list[str]
    """
    ref_idx  = {r: i for i, r in enumerate(entry['refs'])}
    axis_idx = {l: i for i, l in enumerate(entry['axis'])}
    data     = entry['data']

    mat = np.full((len(common_refs), len(common_axis)), np.nan)
    for ri, ref in enumerate(common_refs):
        r = ref_idx[ref]
        for ci, lbl in enumerate(common_axis):
            c     = axis_idx[lbl]
            match = float(data[0, r, c])
            mm    = float(data[1, r, c])
            tot   = match + mm
            if tot > 0:
                mat[ri, ci] = mm / tot
    return mat


def _align_no_base(mask, refs, axis, common_refs, common_axis):
    """
    Re-align a no_base_mask to a common_refs × common_axis subset.
    Returns None if mask is None (backward compat for old NPZ files).
    """
    if mask is None:
        return None
    ref_idx  = {r: i for i, r in enumerate(refs)}
    axis_idx = {l: i for i, l in enumerate(axis)}
    out = np.zeros((len(common_refs), len(common_axis)), dtype=bool)
    for ri, ref in enumerate(common_refs):
        r = ref_idx[ref]
        for ci, lbl in enumerate(common_axis):
            c = axis_idx[lbl]
            out[ri, ci] = bool(mask[r, c])
    return out


def _common_refs_and_axis(loaded):
    """Intersect refs and axis across all loaded entries, preserving first-file order."""
    stems      = list(loaded.keys())
    first_axis = loaded[stems[0]]['axis']
    first_refs = loaded[stems[0]]['refs']
    common_axis = [l for l in first_axis if all(l in set(loaded[s]['axis']) for s in stems)]
    common_refs = [r for r in first_refs if all(r in set(loaded[s]['refs']) for s in stems)]
    return common_refs, common_axis


def _plot_delta_pairs(loaded, common_refs, common_axis, base, ext,
                      ylabel, title, show_insertions, dpi):
    """Compute and save one delta PNG per pair in loaded."""
    for stemA, stemB in itertools.combinations(loaded.keys(), 2):
        rate_A = _mismatch_rate(loaded[stemA], common_refs, common_axis)
        rate_B = _mismatch_rate(loaded[stemB], common_refs, common_axis)
        delta  = rate_A - rate_B

        # Dot if EITHER tRNA lacks the position entirely
        mask_A = _align_no_base(
            loaded[stemA].get('no_base'), loaded[stemA]['refs'],
            loaded[stemA]['axis'], common_refs, common_axis)
        mask_B = _align_no_base(
            loaded[stemB].get('no_base'), loaded[stemB]['refs'],
            loaded[stemB]['axis'], common_refs, common_axis)
        if mask_A is not None and mask_B is not None:
            combined_mask = mask_A | mask_B
        else:
            combined_mask = None

        pair_title = title if title else f"Delta: {stemA} \u2212 {stemB}"
        out_path   = f"{base}_{stemA}_vs_{stemB}{ext}"
        _plot_aligned(delta, common_refs, common_axis, out_path,
                      vmin=-1, vmax=1, cmap_override=_DELTA_CMAP,
                      ylabel=ylabel, title=pair_title,
                      show_insertions=show_insertions, dpi=dpi,
                      no_base_mask=combined_mask)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def plot(arrays, sprinzl_axis, ref_to_sprinzl, output_path,
         palette='light-high',
         ylabel='tRNA Reference Names',
         title='tRNA Alignment Pileup Heatmap',
         show_insertions=False,
         dpi=300):
    """
    Generate and save a heatmap from raw pileup arrays mapped to Sprinzl coordinates.

    Parameters
    ----------
    arrays : dict
        Output of pileup_engine.pileup(). Keys are reference names (str),
        values are np.ndarray of shape (4, L):
            row 0 = match, 1 = mismatch, 2 = insertion, 3 = deletion
    sprinzl_axis : list[str]
        Ordered Sprinzl labels for the shared x-axis.
    ref_to_sprinzl : dict[str, list[str]]
        {ref_name: [sprinzl_label per ref_pos]} from get_sprinzl_mapping().
    output_path : str
        Destination file. Extension determines format (.png, .pdf, .svg).
    palette, ylabel, title, show_insertions : see module docstring.
    """
    ref_names = sorted(arrays.keys())
    n_cols    = len(sprinzl_axis)
    sprinzl_index = {label: i for i, label in enumerate(sprinzl_axis)}

    # TODO: swap in a different metric here if needed.
    # Current metric: mismatch rate = mismatch / (match + mismatch).
    matrix       = np.full((len(ref_names), n_cols), np.nan)
    no_base_mask = np.ones((len(ref_names), n_cols), dtype=bool)

    for row_i, name in enumerate(ref_names):
        arr      = arrays[name]
        match    = arr[0].astype(float)
        mismatch = arr[1].astype(float)
        total    = match + mismatch

        for ref_pos, label in enumerate(ref_to_sprinzl.get(name, [])):
            col = sprinzl_index.get(label)
            if col is None:
                continue
            no_base_mask[row_i, col] = False
            if total[ref_pos] > 0:
                matrix[row_i, col] = mismatch[ref_pos] / total[ref_pos]

    _plot_aligned(matrix, ref_names, sprinzl_axis, output_path,
                  palette=palette, ylabel=ylabel, title=title,
                  show_insertions=show_insertions, dpi=dpi,
                  no_base_mask=no_base_mask)


def plot_from_npz(npz_path, output_path,
                  palette='light-high',
                  ylabel='tRNA Reference Names',
                  title='tRNA Alignment Pileup Heatmap',
                  show_insertions=False,
                  dpi=300):
    """
    Generate and save a heatmap from a pre-computed NPZ file.

    The NPZ must have been produced by save_pileup() and contain:
      data          shape (4, n_refs, n_sprinzl)
      sprinzl_axis  shape (n_sprinzl,)
      ref_names     shape (n_refs,)
      no_base_mask  shape (n_refs, n_sprinzl)  [added in v0.2; optional for old files]

    palette, ylabel, title, show_insertions : see module docstring.
    """
    f            = np.load(npz_path, allow_pickle=True)
    data         = f['data']
    sprinzl_axis = list(f['sprinzl_axis'])
    ref_names    = list(f['ref_names'])
    no_base_mask = f['no_base_mask'] if 'no_base_mask' in f else None

    match    = data[0].astype(float)
    mismatch = data[1].astype(float)
    total    = match + mismatch

    # TODO: swap in a different metric here if needed.
    with np.errstate(invalid='ignore', divide='ignore'):
        matrix = np.where(total > 0, mismatch / total, np.nan)

    _plot_aligned(matrix, ref_names, sprinzl_axis, output_path,
                  palette=palette, ylabel=ylabel, title=title,
                  show_insertions=show_insertions, dpi=dpi,
                  no_base_mask=no_base_mask)


def delta_heatmap(npz_paths, output_prefix,
                  palette='light-high',
                  ylabel='tRNA Reference Names',
                  title=None,
                  show_insertions=False,
                  dpi=300):
    """
    Compute and save pairwise delta heatmaps across a list of NPZ files.

    For each pair (A, B): delta = mismatch_rate_A − mismatch_rate_B.
    NaN where either sample has no coverage. Black dot where either tRNA
    has no base at that Sprinzl position.

    Output files are named:
      {prefix_base}_{stemA}_vs_{stemB}{ext}
    where stem = filename without directory or extension, and ext is taken
    from output_prefix (defaults to '.png').

    Parameters
    ----------
    npz_paths : list[str]
        Paths to pre-computed NPZ files.
    output_prefix : str
        File path prefix (with optional extension, e.g. 'results/delta' or
        'results/delta.pdf').
    palette, ylabel, show_insertions : see module docstring. palette is unused
        for delta plots (always uses Spectral_r diverging colormap).
    title : str or None
        If None, auto-generated as 'Delta: {stemA} − {stemB}' per pair.
    """
    base, ext = os.path.splitext(output_prefix)
    if not ext:
        ext = '.png'

    loaded = {}
    for path in npz_paths:
        f    = np.load(path, allow_pickle=True)
        stem = os.path.splitext(os.path.basename(path))[0]
        loaded[stem] = {
            'data':    f['data'],
            'axis':    list(f['sprinzl_axis']),
            'refs':    list(f['ref_names']),
            'no_base': f['no_base_mask'] if 'no_base_mask' in f else None,
        }

    common_refs, common_axis = _common_refs_and_axis(loaded)

    if not common_axis:
        raise ValueError("No common Sprinzl positions found across the provided NPZ files.")
    if not common_refs:
        raise ValueError("No common reference sequences found across the provided NPZ files.")

    _plot_delta_pairs(loaded, common_refs, common_axis, base, ext,
                      ylabel, title, show_insertions, dpi)


def delta_from_arrays(arrays_by_stem, sprinzl_axis, ref_to_sprinzl,
                      output_prefix,
                      palette='light-high',
                      ylabel='tRNA Reference Names',
                      title=None,
                      show_insertions=False,
                      dpi=300):
    """
    Compute pairwise delta heatmaps directly from pileup arrays (no NPZ required).

    Equivalent to delta_heatmap() but accepts raw pileup data instead of NPZ
    paths — used by the run-delta subcommand to skip intermediate disk I/O.

    Parameters
    ----------
    arrays_by_stem : dict[str, dict[str, np.ndarray]]
        {stem: arrays} where arrays is the dict returned by pileup_engine.pileup()
        and stem is the label used in output filenames and plot titles.
    sprinzl_axis : list[str]
        Shared Sprinzl x-axis from get_sprinzl_mapping().
    ref_to_sprinzl : dict[str, list[str]]
        Per-reference Sprinzl mapping from get_sprinzl_mapping().
    output_prefix : str
        File path prefix (with optional extension). Produces
        {prefix_base}_{stemA}_vs_{stemB}{ext} per pair.
    palette, ylabel, title, show_insertions : see module docstring.
    """
    base, ext = os.path.splitext(output_prefix)
    if not ext:
        ext = '.png'

    loaded = {}
    for stem, arrays in arrays_by_stem.items():
        data, ref_names, no_base_mask = _build_aligned(arrays, sprinzl_axis, ref_to_sprinzl)
        loaded[stem] = {
            'data':    data,
            'axis':    list(sprinzl_axis),
            'refs':    list(ref_names),
            'no_base': no_base_mask,
        }

    common_refs, common_axis = _common_refs_and_axis(loaded)

    if not common_axis:
        raise ValueError("No common Sprinzl positions found across the provided BAM files.")
    if not common_refs:
        raise ValueError("No common reference sequences found across the provided BAM files.")

    _plot_delta_pairs(loaded, common_refs, common_axis, base, ext,
                      ylabel, title, show_insertions, dpi)


def plot_sprinzl_coverage(sprinzl_axis, ref_to_sprinzl, output_path,
                          palette='light-high',
                          ylabel='tRNA Reference Names',
                          title='tRNA Sprinzl Coverage',
                          show_insertions=False,
                          dpi=300):
    """
    Render a presence/absence heatmap of Sprinzl coordinate coverage.

    No pileup data required. Driven entirely by cmalign output (ref_to_sprinzl).

    Colored cell  : this tRNA has a base at this Sprinzl position.
    Black dot (●) : cmalign assigned no base here for this tRNA.

    Parameters
    ----------
    sprinzl_axis : list[str]
        Ordered Sprinzl labels from get_sprinzl_mapping().
    ref_to_sprinzl : dict[str, list[str]]
        Per-reference Sprinzl mapping from get_sprinzl_mapping().
    output_path : str
        Destination file. Extension determines format (.png, .pdf, .svg).
    palette, ylabel, title, show_insertions : see module docstring.
    """
    ref_names     = sorted(ref_to_sprinzl.keys())
    n_cols        = len(sprinzl_axis)
    sprinzl_index = {label: i for i, label in enumerate(sprinzl_axis)}
    matrix        = np.full((len(ref_names), n_cols), np.nan)
    no_base_mask  = np.ones((len(ref_names), n_cols), dtype=bool)

    for row_i, name in enumerate(ref_names):
        for label in ref_to_sprinzl.get(name, []):
            col = sprinzl_index.get(label)
            if col is not None:
                matrix[row_i, col]       = 1.0
                no_base_mask[row_i, col] = False

    _plot_aligned(matrix, ref_names, sprinzl_axis, output_path,
                  palette=palette, ylabel=ylabel, title=title,
                  show_insertions=show_insertions, dpi=dpi,
                  no_base_mask=no_base_mask,
                  cbar_label='Sprinzl position covered')


def save_pileup(arrays, sprinzl_axis, ref_to_sprinzl, base_path):
    """
    Save pileup data mapped to Sprinzl coordinates as both NPZ and TSV.

    Always writes two files:
      {base_path}.npz  — numpy bundle for future plotting with plot_from_npz()
      {base_path}.tsv  — long-format table for human inspection

    Insertion positions (e.g. '36i1') are included in both files regardless
    of the --include-insertions flag; that flag only affects what is rendered
    in the heatmap plot.

    TSV schema (index = seqname + sprinzl_position):
      match, mismatch, deletion, insertion,
      accuracy      = match / (match + mismatch + deletion + insertion),
      mismatch_rate = mismatch / (match + mismatch)

    NPZ schema:
      data          shape (4, n_refs, n_sprinzl)
                    axis 0 = [match, mismatch, insertion, deletion]
      sprinzl_axis  shape (n_sprinzl,)
      ref_names     shape (n_refs,)
      no_base_mask  shape (n_refs, n_sprinzl) — True = no base at that position
    """
    aligned, ref_names, no_base_mask = _build_aligned(arrays, sprinzl_axis, ref_to_sprinzl)

    # --- NPZ ---
    npz_path = base_path + '.npz'
    np.savez(
        npz_path,
        data=aligned,
        sprinzl_axis=np.array(sprinzl_axis),
        ref_names=np.array(ref_names),
        no_base_mask=no_base_mask,
    )
    print(f"Pileup NPZ saved to {npz_path} (data shape: {aligned.shape})")

    # --- TSV ---
    tsv_path = base_path + '.tsv'
    rows = []
    for row_i, name in enumerate(ref_names):
        for col_i, label in enumerate(sprinzl_axis):
            match     = aligned[0, row_i, col_i]
            mismatch  = aligned[1, row_i, col_i]
            insertion = aligned[2, row_i, col_i]
            deletion  = aligned[3, row_i, col_i]

            denom_acc = match + mismatch + deletion + insertion
            accuracy  = (match / denom_acc) if denom_acc > 0 else float('nan')

            denom_mm = match + mismatch
            mm_rate  = (mismatch / denom_mm) if denom_mm > 0 else float('nan')

            rows.append({
                'seqname':          name,
                'sprinzl_position': label,
                'match':            int(match),
                'mismatch':         int(mismatch),
                'deletion':         int(deletion),
                'insertion':        int(insertion),
                'accuracy':         accuracy,
                'mismatch_rate':    mm_rate,
            })

    df = pd.DataFrame(rows).set_index(['seqname', 'sprinzl_position'])
    df.to_csv(tsv_path, sep='\t')
    print(f"Pileup TSV saved to {tsv_path}")
