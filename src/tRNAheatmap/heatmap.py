"""
Heatmap generation and pileup data export for tRNA pileup data.

Public API
----------
  plot(sprinzl_rates, sprinzl_axis, ref_names, no_base_sets, output_path, **kwargs)
      Render one heatmap from Sprinzl-keyed mismatch-rate dicts.

  delta(sprinzl_rates_by_condition, sprinzl_axis, output_prefix, **kwargs)
      Compute pairwise delta heatmaps from per-condition Sprinzl-keyed rate dicts.

  plot_sprinzl_coverage(sprinzl_axis, ref_to_sprinzl, output_path, **kwargs)
      Presence/absence heatmap showing which Sprinzl positions each tRNA covers.
      No pileup data required — driven by cmalign output alone.

  save_pileup(arrays, sprinzl_axis, ref_to_sprinzl, base_path)
      Save raw pileup counts (match/mm/ins/del + derived rates) as <base_path>.tsv.
      Includes has_base column for round-trip replotting.

  save_rates(sprinzl_rates, sprinzl_axis, ref_names, no_base_sets, base_path, std_dict=None)
      Save mismatch rates (with optional per-position std dev) as <base_path>.tsv.
      Includes has_base column for round-trip replotting.

  load_tsv(path)
      Read a TSV written by save_pileup or save_rates and return Sprinzl-keyed
      rate dicts ready for plot() or delta().

Canonical rate representation
-----------------------------
  All plotting functions use Sprinzl-keyed rate dicts:
    sprinzl_rates : dict[str, dict[str, float]]
      {ref_name: {sprinzl_label: mismatch_rate}}  NaN = no coverage
    no_base_sets  : dict[str, set[str]]
      {ref_name: {sprinzl_label, ...}}  positions where tRNA has no base (black dots)
  BAM-derived rates are projected to this form by pipeline.project_to_sprinzl()
  before reaching the heatmap functions.

Shared keyword arguments
------------------------
  palette          : str   'light-high' (default) | 'dark-high' | 'viridis'
  ylabel           : str   Y-axis label. Default 'tRNA Reference Names'.
  title            : str   Figure title. Default 'tRNA Alignment Pileup Heatmap'.
  show_insertions  : bool  Show insertion columns (Xi1, Xi2 ...) in the plot.
                           They are always saved in TSV. Default False.

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
from mpl_toolkits.axes_grid1 import make_axes_locatable

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
                  seq_annotations=None,
                  seq_fontsize=5.0,
                  mod_map=None,
                  cell_size=0.25,
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
        (e.g. '36i1') are stripped before rendering. They remain in TSV.
    vmin, vmax : float
        Color scale limits. Use vmin=-1, vmax=1 for delta plots.
    cmap_override : str or None
        If set, overrides palette (used for delta plots).
    no_base_mask : np.ndarray (bool), shape (n_refs, n_sprinzl)
        True where the tRNA has no base at that Sprinzl position (cmalign gap).
        Black dots are drawn only at these cells.
    cbar_label : str or None
        Override the colorbar axis label. Defaults to 'Mismatch rate' or
        'Δ Mismatch rate' based on vmin.
    """
    # --- Filter insertion columns ---
    if not show_insertions:
        keep = [i for i, lbl in enumerate(sprinzl_axis)
                if not re.search(r'i\d+', lbl)]
        sprinzl_axis = [sprinzl_axis[i] for i in keep]
        matrix = matrix[:, keep]
        if no_base_mask is not None:
            no_base_mask = no_base_mask[:, keep]
        if seq_annotations is not None:
            seq_annotations = seq_annotations[:, keep]

    n_refs, n_cols = matrix.shape
    cmap_name = cmap_override if cmap_override else PALETTES.get(palette, 'hot')
    cmap_obj = plt.get_cmap(cmap_name).copy()
    cmap_obj.set_bad(alpha=0)   # NaN cells: fully transparent

    # --- Figure sizing: aim for roughly square cells ---
    fig_w = max(6.0, n_cols * cell_size + 3.0)
    fig_h = max(3.0, n_refs * cell_size + 2.0)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # pcolormesh renders vector rectangles in PDF (no raster interpolation blur)
    im = ax.pcolormesh(matrix, cmap=cmap_obj, vmin=vmin, vmax=vmax)
    ax.set_xlim(0, n_cols)
    ax.set_ylim(n_refs, 0)   # top-to-bottom row order, matching imshow convention

    # --- Black dots: no-base positions only ---
    if no_base_mask is not None:
        dot_rows, dot_cols = np.where(no_base_mask)
        if dot_rows.size > 0:
            # pcolormesh cell centres are at col+0.5, row+0.5
            ax.scatter(dot_cols + 0.5, dot_rows + 0.5, s=20, c='black', marker='o',
                       zorder=2, linewidths=0)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size=f"{200 / n_cols:.2f}%", pad="1%")
    cbar = fig.colorbar(im, cax=cax)
    cbar.ax.tick_params()
    if cbar_label is None:
        cbar_label = 'Mismatch rate' if vmin == 0 else 'Δ Mismatch rate'
    cbar.ax.set_ylabel(cbar_label, labelpad=10, rotation=270, va='bottom')

    # pcolormesh cell centres at i+0.5
    ax.set_yticks([i + 0.5 for i in range(n_refs)])
    ax.set_yticklabels(ref_names)
    ax.tick_params(axis='y', length=4)

    ax.set_xticks([i + 0.5 for i in range(n_cols)])
    ax.set_xticklabels(sprinzl_axis, rotation=90)
    ax.tick_params(axis='x', length=6)

    ax.set_xlabel('Sprinzl position')
    ax.set_ylabel(ylabel)
    ax.set_title(title, pad=7)

    # Build a set of (row, col) cells that have a modification symbol so the
    # sequence annotation loop can skip those cells (mods take priority).
    mod_cells = set()
    if mod_map:
        for col_i, label in enumerate(sprinzl_axis):
            for row_i, ref_name in enumerate(ref_names):
                symbol = mod_map.get(ref_name, {}).get(label)
                if not symbol:
                    continue
                if no_base_mask is not None and no_base_mask[row_i, col_i]:
                    continue
                val = matrix[row_i, col_i]
                if np.isnan(val):
                    continue
                norm = np.clip((val - vmin) / (vmax - vmin), 0, 1)
                r, g, b, _ = cmap_obj(norm)
                text_color = 'white' if (0.299*r + 0.587*g + 0.114*b) < 0.5 else 'black'
                ax.text(col_i + 0.5, row_i + 0.5, symbol,
                        ha='center', va='center',
                        fontsize=seq_fontsize, color=text_color,
                        fontweight='bold', zorder=4)
                mod_cells.add((row_i, col_i))

    if seq_annotations is not None:
        for row_i in range(n_refs):
            for col_i in range(n_cols):
                if (row_i, col_i) in mod_cells:
                    continue
                letter = seq_annotations[row_i, col_i]
                if letter:
                    ax.text(col_i + 0.5, row_i + 0.5, letter,
                            ha='center', va='center',
                            fontsize=seq_fontsize, color='white', fontweight='bold',
                            zorder=3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, transparent=True)
    plt.close(fig)
    print(f"Heatmap saved to {output_path}")


def _build_count_matrix(arrays, sprinzl_axis, ref_to_sprinzl):
    """
    Map per-reference 4-row pileup count arrays to Sprinzl coordinate columns.

    Used only by save_pileup() — plotting/delta operate on rate matrices.

    Returns
    -------
    aligned : np.ndarray, shape (4, n_refs, n_sprinzl)
        axis 0 = [match, mismatch, insertion, deletion]
    ref_names : list[str]
    """
    ref_names     = list(arrays.keys())
    n_sprinzl     = len(sprinzl_axis)
    sprinzl_index = {label: i for i, label in enumerate(sprinzl_axis)}
    aligned       = np.zeros((4, len(ref_names), n_sprinzl), dtype=np.float64)

    for row_i, name in enumerate(ref_names):
        arr = arrays[name]
        for ref_pos, label in enumerate(ref_to_sprinzl.get(name, [])):
            col = sprinzl_index.get(label)
            if col is not None and ref_pos < arr.shape[1]:
                aligned[:, row_i, col] = arr[:, ref_pos]

    return aligned, ref_names


def _build_rate_matrix(rate_dict, sprinzl_axis, ref_to_sprinzl):
    """
    Map per-reference 1D mismatch-rate arrays to Sprinzl coordinate columns.

    Parameters
    ----------
    rate_dict : dict[str, np.ndarray]
        {ref_name: 1D float64 array of mismatch rates, NaN where uncovered}.
    sprinzl_axis : list[str]
    ref_to_sprinzl : dict[str, list[str]]

    Returns
    -------
    matrix : np.ndarray, shape (n_refs, n_sprinzl)
        Mismatch rates in Sprinzl space. NaN = no coverage.
    ref_names : list[str]
    no_base_mask : np.ndarray (bool), shape (n_refs, n_sprinzl)
        True where the tRNA has no base at that Sprinzl position (cmalign gap).
    """
    ref_names     = list(rate_dict.keys())
    n_sprinzl     = len(sprinzl_axis)
    sprinzl_index = {label: i for i, label in enumerate(sprinzl_axis)}
    matrix        = np.full((len(ref_names), n_sprinzl), np.nan)
    no_base_mask  = np.ones((len(ref_names), n_sprinzl), dtype=bool)

    for row_i, name in enumerate(ref_names):
        rates = rate_dict[name]
        for ref_pos, label in enumerate(ref_to_sprinzl.get(name, [])):
            col = sprinzl_index.get(label)
            if col is None:
                continue
            no_base_mask[row_i, col] = False
            if ref_pos < len(rates) and not np.isnan(rates[ref_pos]):
                matrix[row_i, col] = rates[ref_pos]

    return matrix, ref_names, no_base_mask


def _build_matrix_from_sprinzl_rates(sprinzl_rates, sprinzl_axis, ref_names, no_base_sets):
    """
    Build a plot matrix directly from Sprinzl-keyed rate dicts.

    Parameters
    ----------
    sprinzl_rates : dict[str, dict[str, float]]
        {ref_name: {sprinzl_label: mismatch_rate}}. NaN = no coverage.
    sprinzl_axis  : list[str]
    ref_names     : list[str]  — ordered row list (preserves include-refs order)
    no_base_sets  : dict[str, set[str]]
        {ref_name: {sprinzl_label, ...}} — positions with no base (black dots).

    Returns
    -------
    matrix       : np.ndarray, shape (n_refs, n_sprinzl)
    no_base_mask : np.ndarray (bool), shape (n_refs, n_sprinzl)
    """
    sprinzl_index = {lbl: i for i, lbl in enumerate(sprinzl_axis)}
    n_refs = len(ref_names)
    n_cols = len(sprinzl_axis)
    matrix       = np.full((n_refs, n_cols), np.nan)
    no_base_mask = np.zeros((n_refs, n_cols), dtype=bool)

    for row_i, name in enumerate(ref_names):
        ref_rates = sprinzl_rates.get(name, {})
        no_base   = no_base_sets.get(name, set())
        for lbl, col in sprinzl_index.items():
            if lbl in no_base:
                no_base_mask[row_i, col] = True
            elif lbl in ref_rates:
                val = ref_rates[lbl]
                if not np.isnan(val):
                    matrix[row_i, col] = val
    return matrix, no_base_mask


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def plot(sprinzl_rates, sprinzl_axis, ref_names, no_base_sets, output_path,
         palette='light-high',
         ylabel='tRNA Reference Names',
         title='tRNA Alignment Pileup Heatmap',
         show_insertions=False,
         mod_map=None,
         cell_size=0.25,
         dpi=300):
    """
    Render a single heatmap from Sprinzl-keyed mismatch-rate dicts.

    Parameters
    ----------
    sprinzl_rates : dict[str, dict[str, float]]
        {ref_name: {sprinzl_label: mismatch_rate}}. NaN = no coverage.
        Use pipeline.project_to_sprinzl() for BAM-derived rates, or
        load_tsv() for TSV-sourced rates.
    sprinzl_axis : list[str]
        Ordered Sprinzl labels for the shared x-axis.
    ref_names : list[str]
        Ordered list of reference names for the y-axis (preserves include-refs order).
    no_base_sets : dict[str, set[str]]
        {ref_name: {sprinzl_label, ...}} — positions shown as black dots.
    output_path : str
        Destination file. Extension determines format (.png, .pdf, .svg).
    palette, ylabel, title, show_insertions : see module docstring.
    """
    matrix, no_base_mask = _build_matrix_from_sprinzl_rates(
        sprinzl_rates, sprinzl_axis, ref_names, no_base_sets)
    _plot_aligned(matrix, ref_names, sprinzl_axis, output_path,
                  palette=palette, ylabel=ylabel, title=title,
                  show_insertions=show_insertions, dpi=dpi,
                  cell_size=cell_size, mod_map=mod_map, no_base_mask=no_base_mask)


def delta(sprinzl_rates_by_condition, sprinzl_axis, output_prefix,
          no_base_sets_by_condition=None,
          palette='light-high',
          ylabel='tRNA Reference Names',
          title=None,
          show_insertions=False,
          mod_map=None,
          cell_size=0.25,
          dpi=300):
    """
    Compute pairwise delta heatmaps from per-condition Sprinzl-keyed rate dicts.

    Each pair (condA, condB) produces one file named
    {output_prefix}_{condA}_vs_{condB}{ext}. The reference set is intersected
    across conditions; only refs present in all conditions are compared.

    Parameters
    ----------
    sprinzl_rates_by_condition : dict[str, dict[str, dict[str, float]]]
        {condition_name: {ref_name: {sprinzl_label: mismatch_rate}}}.
    sprinzl_axis : list[str]
        Shared Sprinzl x-axis.
    output_prefix : str
        File path prefix (with optional extension).
    no_base_sets_by_condition : dict[str, dict[str, set[str]]] or None
        {condition_name: {ref_name: {sprinzl_label, ...}}}. Black dot appears
        if EITHER condition has no base at that position.
    title : str or None
        If None, defaults to "Delta: {condA} − {condB}" per pair.
    palette, ylabel, show_insertions, mod_map, cell_size, dpi : see plot().
    """
    base, ext = os.path.splitext(output_prefix)
    if not ext:
        ext = '.pdf'
    if no_base_sets_by_condition is None:
        no_base_sets_by_condition = {}

    cond_names = list(sprinzl_rates_by_condition.keys())
    # Intersect ref sets across all conditions, preserving order of first condition
    first_refs = list(sprinzl_rates_by_condition[cond_names[0]].keys())
    common_set = set.intersection(*(set(sprinzl_rates_by_condition[c].keys())
                                    for c in cond_names))
    common_refs = [r for r in first_refs if r in common_set]
    if not common_refs:
        raise ValueError("No common reference sequences found across conditions.")

    for condA, condB in itertools.combinations(cond_names, 2):
        mA, maskA = _build_matrix_from_sprinzl_rates(
            sprinzl_rates_by_condition[condA], sprinzl_axis, common_refs,
            no_base_sets_by_condition.get(condA, {}))
        mB, maskB = _build_matrix_from_sprinzl_rates(
            sprinzl_rates_by_condition[condB], sprinzl_axis, common_refs,
            no_base_sets_by_condition.get(condB, {}))

        delta_mat     = mA - mB
        combined_mask = maskA | maskB

        pair_title = title if title else f"Delta: {condA} − {condB}"
        out_path   = f"{base}_{condA}_vs_{condB}{ext}"
        _plot_aligned(delta_mat, common_refs, sprinzl_axis, out_path,
                      vmin=-1, vmax=1, cmap_override=_DELTA_CMAP,
                      ylabel=ylabel, title=pair_title,
                      show_insertions=show_insertions, dpi=dpi,
                      cell_size=cell_size, mod_map=mod_map,
                      no_base_mask=combined_mask)


def plot_sprinzl_coverage(sprinzl_axis, ref_to_sprinzl, output_path,
                          palette='light-high',
                          ylabel='tRNA Reference Names',
                          title='tRNA Sprinzl Coverage',
                          show_insertions=False,
                          ref_to_seq=None,
                          seq_fontsize=5.0,
                          mod_map=None,
                          cell_size=0.25,
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
    ref_names     = list(ref_to_sprinzl.keys())
    n_cols        = len(sprinzl_axis)
    sprinzl_index = {label: i for i, label in enumerate(sprinzl_axis)}
    matrix        = np.full((len(ref_names), n_cols), np.nan)
    no_base_mask  = np.ones((len(ref_names), n_cols), dtype=bool)

    seq_annotations = None
    if ref_to_seq is not None:
        seq_annotations = np.full((len(ref_names), n_cols), '', dtype=object)

    for row_i, name in enumerate(ref_names):
        labels = ref_to_sprinzl.get(name, [])
        seq    = ref_to_seq.get(name, '') if ref_to_seq is not None else ''
        for seq_pos, label in enumerate(labels):
            col = sprinzl_index.get(label)
            if col is not None:
                matrix[row_i, col]       = 1.0
                no_base_mask[row_i, col] = False
                if seq_annotations is not None and seq_pos < len(seq):
                    seq_annotations[row_i, col] = seq[seq_pos].upper()

    _plot_aligned(matrix, ref_names, sprinzl_axis, output_path,
                  palette=palette, ylabel=ylabel, title=title,
                  show_insertions=show_insertions, dpi=dpi,
                  cell_size=cell_size, mod_map=mod_map, no_base_mask=no_base_mask,
                  seq_annotations=seq_annotations, seq_fontsize=seq_fontsize,
                  cbar_label='Sprinzl position covered')


def save_pileup(arrays, sprinzl_axis, ref_to_sprinzl, base_path):
    """
    Save raw pileup counts mapped to Sprinzl coordinates as {base_path}.tsv.

    Insertion positions (e.g. '36i1') are included regardless of the
    --include-insertions flag; that flag only affects what is rendered in the plot.

    TSV schema (index = seqname + sprinzl_position):
      match, mismatch, deletion, insertion,
      accuracy      = match / (match + mismatch + deletion + insertion),
      mismatch_rate = mismatch / (match + mismatch + deletion + insertion)
    """
    aligned, ref_names = _build_count_matrix(arrays, sprinzl_axis, ref_to_sprinzl)
    ref_label_sets = {name: set(ref_to_sprinzl.get(name, [])) for name in ref_names}

    rows = []
    for row_i, name in enumerate(ref_names):
        for col_i, label in enumerate(sprinzl_axis):
            match     = aligned[0, row_i, col_i]
            mismatch  = aligned[1, row_i, col_i]
            insertion = aligned[2, row_i, col_i]
            deletion  = aligned[3, row_i, col_i]

            denom_acc = match + mismatch + deletion + insertion
            accuracy  = (match / denom_acc) if denom_acc > 0 else float('nan')

            denom_mm = match + mismatch + insertion + deletion
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
                'has_base':         label in ref_label_sets[name],
            })

    tsv_path = base_path + '.tsv'
    df = pd.DataFrame(rows).set_index(['seqname', 'sprinzl_position'])
    df.to_csv(tsv_path, sep='\t')
    print(f"Pileup TSV saved to {tsv_path}")


def save_rates(sprinzl_rates, sprinzl_axis, ref_names, no_base_sets, base_path,
               std_dict=None):
    """
    Save per-position mismatch rates (and optional std deviation) to {base_path}.tsv.

    TSV schema (index = seqname + sprinzl_position):
      mismatch_rate  — mismatch rate (0–1), NaN = no coverage
      std_dev        — std deviation across BAMs (only when std_dict provided)
      has_base       — True if the tRNA has a base at this position; False = cmalign gap
    """
    rows = []
    for name in ref_names:
        ref_rates = sprinzl_rates.get(name, {})
        ref_stds  = (std_dict or {}).get(name, {})
        no_base   = no_base_sets.get(name, set())
        for label in sprinzl_axis:
            row = {
                'seqname':          name,
                'sprinzl_position': label,
                'mismatch_rate':    ref_rates.get(label, float('nan')),
                'has_base':         label not in no_base,
            }
            if std_dict is not None:
                row['std_dev'] = ref_stds.get(label, float('nan'))
            rows.append(row)

    tsv_path = base_path + '.tsv'
    df = pd.DataFrame(rows).set_index(['seqname', 'sprinzl_position'])
    df.to_csv(tsv_path, sep='\t')
    print(f"Rate TSV saved to {tsv_path}")


def load_tsv(path):
    """
    Read a TSV written by save_pileup or save_rates and return Sprinzl-keyed
    rate dicts ready for plot() or delta().

    Detects pileup format (has 'match' column) vs rates format automatically.
    Reads the optional 'has_base' column to reconstruct black-dot positions;
    older TSVs without the column produce empty no_base_sets (no black dots).

    Parameters
    ----------
    path : str

    Returns
    -------
    sprinzl_rates : dict[str, dict[str, float]]
        {ref_name: {sprinzl_label: mismatch_rate}}. NaN = no coverage.
    sprinzl_axis  : list[str]
        Canonical-sorted union of all Sprinzl labels in the file.
    no_base_sets  : dict[str, set[str]]
        {ref_name: {sprinzl_label, ...}} — cmalign-gap positions (black dots).
        Empty per ref when 'has_base' column is absent.
    ref_names     : list[str]
        References in the order they appear in the file (preserves include-refs order).
    """
    from .calculate_tRNA_positions import _sprinzl_sort_key

    df = pd.read_csv(path, sep='\t')
    is_pileup    = 'match' in df.columns
    has_base_col = 'has_base' in df.columns

    sprinzl_rates = {}
    no_base_sets  = {}
    ref_names     = list(dict.fromkeys(df['seqname']))  # preserves order, deduplicates

    for name in ref_names:
        sub = df[df['seqname'] == name]
        rates   = {}
        no_base = set()
        for _, row in sub.iterrows():
            label = row['sprinzl_position']
            if is_pileup:
                total = (row['match'] + row['mismatch']
                         + row['insertion'] + row['deletion'])
                rate = float(row['mismatch'] / total) if total > 0 else float('nan')
            else:
                rate = float(row['mismatch_rate'])
            if has_base_col and not row['has_base']:
                no_base.add(label)
            else:
                rates[label] = rate
        sprinzl_rates[name] = rates
        if no_base:
            no_base_sets[name] = no_base

    all_labels   = {lbl for d in sprinzl_rates.values() for lbl in d}
    all_labels  |= {lbl for s in no_base_sets.values() for lbl in s}
    sprinzl_axis = sorted(all_labels, key=_sprinzl_sort_key)

    return sprinzl_rates, sprinzl_axis, no_base_sets, ref_names
