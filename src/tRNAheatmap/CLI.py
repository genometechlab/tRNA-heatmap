#!/usr/bin/env python3
"""
CLI entry point for tRNA pileup heatmap generation.

Subcommands
-----------
  run      Full pipeline: pileup → Sprinzl mapping → heatmap. With one input
           condition, emits a single heatmap; with ≥2, emits pairwise delta
           heatmaps (and per-condition heatmaps with --individual).
  inspect  Dry-run: visualize Sprinzl coordinate coverage without any pileup.

Examples
--------
  # Single heatmap from one BAM (one --condition group with one BAM)
  python -m tRNAheatmap.CLI run \\
      --ref ref.fa --organism eukaryotic \\
      --condition sample aligned.bam \\
      --output heatmap.pdf --save-df pileup

  # Pairwise delta from two single-BAM conditions
  python -m tRNAheatmap.CLI run \\
      --ref ref.fa --organism eukaryotic \\
      --condition condA condA.bam --condition condB condB.bam \\
      --output results/delta

  # Replicate-aware deltas with --merge-mode equal
  python -m tRNAheatmap.CLI run \\
      --ref ref.fa --organism eukaryotic --merge-mode equal \\
      --condition WT wt1.bam wt2.bam wt3.bam \\
      --condition KO ko1.bam ko2.bam ko3.bam \\
      --output results/delta

  # Inspect Sprinzl coordinate coverage (no BAM or pileup needed)
  python -m tRNAheatmap.CLI inspect --ref ref.fa --organism eukaryotic --output coverage.pdf
"""

import argparse
import os

from .calculate_tRNA_positions import (get_sprinzl_mapping, save_sprinzl_mapping,
                                        load_sprinzl_mapping, build_axis_from_mapping,
                                        _sprinzl_sort_key)
import matplotlib.pyplot as plt
from . import heatmap, pipeline

_PKG_DIR = os.path.dirname(os.path.abspath(__file__))
_CM_DIR  = os.path.join(_PKG_DIR, 'alignment_cm')
BUNDLED_MODELS = {
    'eukaryotic': os.path.join(_CM_DIR, 'euk-num.cm'),
    'prokaryotic': os.path.join(_CM_DIR, 'bact-num.cm'),
    'archaeal':    os.path.join(_CM_DIR, 'arch-num.cm'),
}

PALETTE_CHOICES = list(heatmap.PALETTES.keys())


# ---------------------------------------------------------------------------
# Argparse helpers
# ---------------------------------------------------------------------------

def _add_plot_flags(p):
    """Add shared visualisation flags to a subparser."""
    p.add_argument(
        "--palette",
        choices=PALETTE_CHOICES,
        default='light-high',
        help="Color palette. 'light-high': light = more mismatch (default). "
             "'dark-high': dark = more mismatch. 'viridis': perceptually uniform."
    )
    p.add_argument(
        "--ylabel",
        default="tRNA Reference Names",
        help="Y-axis label. Default: 'tRNA Reference Names'."
    )
    p.add_argument(
        "--title",
        default="tRNA Alignment Pileup Heatmap",
        help="Figure title. Default: 'tRNA Alignment Pileup Heatmap'."
    )
    p.add_argument(
        "--include-insertions",
        action="store_true",
        dest="include_insertions",
        help="Show insertion positions (e.g. 36i1) in the heatmap. "
             "They are always stored in TSV files."
    )
    p.add_argument(
        "--reference-match",
        action="store_true",
        dest="reference_match",
        help="Plot reference match proportion (match / (match+mismatch+ins+del)) "
             "instead of mismatch rate. Colorbar label updates accordingly. "
             "Delta plots become Δ Reference match rate. Note: the TSV "
             "'mismatch_rate' column name is unchanged; with this flag its "
             "values hold the match rate instead."
    )
    p.add_argument(
        "--annotate",
        action="store_true",
        dest="annotate",
        help="Print each cell's numeric value on the heatmap. Delta cells are "
             "shown signed (e.g. +0.42); rate cells unsigned (e.g. 0.42). "
             "Cells with a modification symbol keep the symbol; no-coverage and "
             "no-base cells stay blank."
    )
    p.add_argument(
        "--annotate-fontsize",
        dest="annotate_fontsize",
        type=float,
        default=5.0,
        metavar="FLOAT",
        help="Font size for numbers drawn with --annotate. Default: 5.0."
    )
    p.add_argument(
        "--annotate-threshold",
        dest="annotate_threshold",
        type=float,
        default=0.0,
        metavar="FLOAT",
        help="With --annotate, only print a cell's value when its absolute value "
             "is greater than or equal to this threshold. Default: 0.0 (print all "
             "covered cells, including exact-zero deltas)."
    )
    p.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Output image resolution in DPI. Default: 300."
    )
    p.add_argument(
        "--cell-size",
        dest="cell_size",
        type=float,
        default=0.25,
        metavar="FLOAT",
        help="Inches per heatmap cell for auto figure sizing. Default: 0.25."
    )
    p.add_argument(
        "--style",
        dest="mpl_style",
        default=None,
        metavar="FILE",
        help="Matplotlib style sheet (.mplstyle) to layer on top of the bundled default. "
             "Copy src/tRNAheatmap/plotting_styles/default.mplstyle as a starting point."
    )
    p.add_argument(
        "--outdir", "-O",
        default=None,
        dest="outdir",
        help="Directory for all output files. Created if it does not exist. "
             "Default: current working directory."
    )


def _add_adapter_flags(p):
    """Add mutually-exclusive adapter trimming flags to a subparser."""
    adapt_group = p.add_mutually_exclusive_group()
    adapt_group.add_argument(
        "--detect-adapters",
        action="store_true",
        dest="detect_adapters",
        help="Auto-detect common 5'/3' adapter sequences shared across all tRNAs "
             "in the reference FASTA and trim them before Sprinzl assignment. "
             "Requires ≥2 sequences. Mutually exclusive with --trim-5/--trim-3."
    )
    adapt_group.add_argument(
        "--trim-5",
        type=int,
        default=0,
        dest="trim_5",
        metavar="N",
        help="Trim N bases from the 5' end of each reference sequence before "
             "Sprinzl coordinate assignment. Mutually exclusive with --detect-adapters."
    )
    # --trim-3 is independent of --detect-adapters but grouped with --trim-5
    # Use a separate argument (outside the group) so both 5' and 3' can be set together
    p.add_argument(
        "--trim-3",
        type=int,
        default=0,
        dest="trim_3",
        metavar="N",
        help="Trim N bases from the 3' end of each reference sequence before "
             "Sprinzl coordinate assignment."
    )


def _add_sprinzl_override_flags(p):
    """Add --sprinzl-map and --save-mapping flags to a subparser."""
    p.add_argument(
        "--sprinzl-map",
        dest="sprinzl_map",
        default=None,
        metavar="TSV",
        help="Pre-computed Sprinzl mapping TSV (from --save-mapping). "
             "When combined with --organism/--cm, patches only the listed refs "
             "and leaves others as computed by cmalign. Without --organism/--cm, "
             "skips cmalign entirely (requires a #axis: header and all refs listed)."
    )
    p.add_argument(
        "--save-mapping",
        dest="save_mapping",
        default=None,
        metavar="BASE",
        help="Save the computed Sprinzl mapping to BASE.tsv for inspection or "
             "reuse with --sprinzl-map."
    )


def _add_condition_flags(p):
    """Add --condition, --condition-ref, and --merge-mode flags to a subparser."""
    p.add_argument(
        "--condition", "-C",
        nargs='+',
        action='append',
        dest='condition',
        default=None,
        metavar='NAME_OR_BAM',
        help="Condition group: first token is the condition name, "
             "remaining tokens are sorted BAM paths. "
             "Repeat for multiple conditions: "
             "--condition WT rep1.bam rep2.bam --condition KO rep3.bam rep4.bam. "
             "At least one --condition is required."
    )
    p.add_argument(
        "--condition-ref",
        nargs=2,
        action='append',
        dest='condition_ref',
        default=None,
        metavar=('NAME', 'FASTA'),
        help="Override the reference FASTA for a named condition. NAME must match "
             "one of the --condition names. Repeat for multiple conditions: "
             "--condition-ref WT wt_ref.fa --condition-ref KO ko_ref.fa. "
             "Conditions without a --condition-ref fall back to --ref. "
             "--ref may be omitted entirely if all BAM conditions have an explicit "
             "--condition-ref."
    )
    p.add_argument(
        "--merge-mode",
        choices=['total', 'equal'],
        default='total',
        dest='merge_mode',
        help="How to merge multiple BAMs within a condition. "
             "'total' (default): sum raw counts — higher-depth replicates "
             "contribute proportionally more. "
             "'equal': average per-BAM mismatch rates — each replicate "
             "contributes equally regardless of sequencing depth."
    )


def _add_ref_filter_flags(p):
    """Add mutually-exclusive reference filtering flags to a subparser."""
    ref_group = p.add_mutually_exclusive_group()
    ref_group.add_argument(
        "--include-refs",
        dest="include_refs",
        default=None,
        metavar="FILE",
        help="Plain-text file (one reference name per line). Only the listed "
             "references will appear in the output; all others are dropped. "
             "Mutually exclusive with --exclude-refs."
    )
    ref_group.add_argument(
        "--exclude-refs",
        dest="exclude_refs",
        default=None,
        metavar="FILE",
        help="Plain-text file (one reference name per line). Listed references "
             "are removed from the output; all others are kept. "
             "Mutually exclusive with --include-refs."
    )


def _filter_refs(arrays, args):
    """Thin shim: pull include/exclude paths off args and dispatch to pipeline."""
    return pipeline.apply_ref_filter(
        arrays,
        include_file=getattr(args, 'include_refs', None),
        exclude_file=getattr(args, 'exclude_refs', None),
    )


def _resolve_sprinzl_mapping(args, ref_for_cmalign, parser):
    """
    Compute or load the Sprinzl mapping based on CLI flags.

    Modes
    -----
    --organism/--cm only        : run cmalign, return result directly.
    --sprinzl-map only          : load TSV; requires #axis: header and all refs.
    both                        : run cmalign, then patch with TSV entries.
    neither                     : error.
    """
    has_cm = bool(getattr(args, 'organism', None) or getattr(args, 'cm', None))
    has_map = bool(getattr(args, 'sprinzl_map', None))

    if not has_cm and not has_map:
        parser.error("Provide --organism / --cm, or --sprinzl-map (or both).")

    if has_cm:
        cm_path = BUNDLED_MODELS[args.organism] if args.organism else args.cm
        print(f"Aligning {ref_for_cmalign} to {cm_path} to compute Sprinzl coordinates...")
        sprinzl_axis, ref_to_sprinzl = get_sprinzl_mapping(ref_for_cmalign, cm_path)
    else:
        sprinzl_axis, ref_to_sprinzl = None, {}

    mod_map = {}
    if has_map:
        tsv_mapping, mod_map = load_sprinzl_mapping(args.sprinzl_map)
        if not has_cm:
            ref_to_sprinzl = tsv_mapping
            sprinzl_axis = build_axis_from_mapping(ref_to_sprinzl)
        else:
            # Patch mode: TSV entries overwrite cmalign results for listed refs;
            # rebuild axis so any new labels (e.g. -1) from the TSV are included.
            ref_to_sprinzl.update(tsv_mapping)
            sprinzl_axis = build_axis_from_mapping(ref_to_sprinzl)

    return sprinzl_axis, ref_to_sprinzl, mod_map


def build_parser():
    parser = argparse.ArgumentParser(
        description="tRNA pileup heatmap tool.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest='command', required=True)

    # ------------------------------------------------------------------
    # 'run' subcommand — full pipeline (single heatmap or pairwise deltas)
    # ------------------------------------------------------------------
    run_p = subparsers.add_parser(
        'run',
        help='Full pipeline: pileup + Sprinzl mapping + heatmap. With ≥2 '
             '--condition groups, emits pairwise delta heatmaps.',
    )
    run_p.add_argument(
        "--ref", "-r",
        default=None,
        help="Path to reference tRNA FASTA file. Required when any --condition "
             "specifies BAM files; optional when all conditions are TSV files."
    )
    cm_group = run_p.add_mutually_exclusive_group(required=False)
    cm_group.add_argument(
        "--organism", "-g",
        choices=["eukaryotic", "prokaryotic", "archaeal"],
        dest="organism",
        help="Use a bundled covariance model: eukaryotic (euk-num.cm), "
             "prokaryotic (bact-num.cm), or archaeal (arch-num.cm)."
    )
    cm_group.add_argument(
        "--cm", "-c",
        dest="cm",
        metavar="FILE",
        help="Path to a custom Infernal covariance model (.cm file)."
    )
    run_p.add_argument(
        "--output", "-o",
        default="heatmap.pdf",
        help="With one --condition: output heatmap file (extension determines format). "
             "With ≥2 --condition groups: prefix; pairs become {prefix}_{A}_vs_{B}.{ext} "
             "(default extension .pdf). Default: heatmap.pdf"
    )
    run_p.add_argument(
        "--save-df", "-s",
        default=None,
        dest="save_df",
        help="Base name (no extension) to save pileup/rate data as .tsv file(s). "
             "With one --condition: BASE.tsv. With ≥ 2 --condition groups: "
             "BASE_{condition}.tsv per condition, plus BASE_{condA}_vs_{condB}.tsv "
             "with the pairwise delta values for each pair."
    )
    run_p.add_argument(
        "--individual",
        action="store_true",
        dest="individual",
        help="With ≥2 --condition groups, also emit a per-condition heatmap "
             "({prefix}_{condition}.{ext}) alongside the pairwise deltas."
    )
    run_p.add_argument(
        "--threads", "-t",
        type=int,
        default=1,
        help="Threads for pileup engine. Default: 1"
    )
    run_p.add_argument(
        "--min-q",
        type=int,
        default=0,
        dest="min_q",
        metavar="INT",
        help="Minimum MAPQ for a read to be included in the pileup. "
             "Reads with mapping_quality < min_q are skipped. Default: 0 (no filter)."
    )
    _add_adapter_flags(run_p)
    _add_sprinzl_override_flags(run_p)
    _add_condition_flags(run_p)
    _add_ref_filter_flags(run_p)
    _add_plot_flags(run_p)

    # ------------------------------------------------------------------
    # 'inspect' subcommand — Sprinzl coverage without pileup
    # ------------------------------------------------------------------
    inspect_p = subparsers.add_parser(
        'inspect',
        help='Dry-run: visualize Sprinzl coordinate coverage without any pileup.',
    )
    inspect_p.add_argument(
        "--ref", "-r",
        required=True,
        help="Path to reference tRNA FASTA file."
    )
    cm_group_in = inspect_p.add_mutually_exclusive_group(required=False)
    cm_group_in.add_argument(
        "--organism", "-g",
        choices=["eukaryotic", "prokaryotic", "archaeal"],
        dest="organism",
        help="Use a bundled covariance model: eukaryotic (euk-num.cm), "
             "prokaryotic (bact-num.cm), or archaeal (arch-num.cm)."
    )
    cm_group_in.add_argument(
        "--cm", "-c",
        dest="cm",
        metavar="FILE",
        help="Path to a custom Infernal covariance model (.cm file)."
    )
    inspect_p.add_argument(
        "--output", "-o",
        default="sprinzl_coverage.pdf",
        help="Output heatmap file. Format from extension (.png, .pdf, .svg). "
             "Default: sprinzl_coverage.pdf"
    )
    inspect_p.add_argument(
        "--include-sequence",
        dest="include_sequence",
        action="store_true",
        default=False,
        help="Write the nucleotide letter in white on each covered cell.",
    )
    inspect_p.add_argument(
        "--seq-fontsize",
        dest="seq_fontsize",
        type=float,
        default=5.0,
        metavar="FLOAT",
        help="Font size for nucleotide letters drawn with --include-sequence. Default: 5.0.",
    )
    _add_ref_filter_flags(inspect_p)
    _add_adapter_flags(inspect_p)
    _add_sprinzl_override_flags(inspect_p)
    _add_plot_flags(inspect_p)

    return parser


def _resolve_output(path, outdir):
    """Redirect path into outdir, creating the directory if needed."""
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        return os.path.join(outdir, os.path.basename(path))
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    return path


def _trim_context_for_ref(ref_path, args, parser):
    """
    Return a pipeline.adapter_trimmed_ref context manager for an explicit ref path.

    Adapter flags (--detect-adapters, --trim-5, --trim-3) are read from args,
    but the FASTA path comes from ref_path rather than args.ref so that
    per-condition references work correctly.
    """
    detect = getattr(args, 'detect_adapters', False)
    trim5  = getattr(args, 'trim_5', 0)
    trim3  = getattr(args, 'trim_3', 0)
    if detect and (trim5 or trim3):
        parser.error("--detect-adapters cannot be combined with --trim-5 or --trim-3.")
    if detect:
        seqs = pipeline.read_fasta(ref_path)
        if len(seqs) < 2:
            parser.error(f"--detect-adapters requires ≥2 sequences in {ref_path!r}.")
    return pipeline.adapter_trimmed_ref(
        ref_path, trim_5=trim5, trim_3=trim3, detect=detect)


def _adapter_trim_context(args, parser):
    """Thin wrapper used by the inspect subcommand (single global --ref)."""
    return _trim_context_for_ref(args.ref, args, parser)


def _build_condition_refs(args, bam_conditions, parser):
    """
    Build {cond_name: fasta_path} for every BAM condition.

    --condition-ref NAME FASTA entries override --ref for specific conditions.
    --ref is used as the fallback for any condition without an explicit override.
    Raises parser.error if any BAM condition has no resolvable reference.
    """
    overrides = {}
    for pair in (args.condition_ref or []):
        name, fasta = pair
        if name not in bam_conditions:
            parser.error(
                f"--condition-ref: '{name}' does not match any BAM --condition. "
                f"Known BAM conditions: {list(bam_conditions)}"
            )
        if name in overrides:
            parser.error(f"--condition-ref: duplicate entry for condition '{name}'.")
        if not os.path.isfile(fasta):
            parser.error(f"--condition-ref: FASTA not found: {fasta!r}")
        overrides[name] = fasta

    condition_refs = {}
    missing = []
    for cond_name in bam_conditions:
        if cond_name in overrides:
            condition_refs[cond_name] = overrides[cond_name]
        elif args.ref:
            condition_refs[cond_name] = args.ref
        else:
            missing.append(cond_name)

    if missing:
        parser.error(
            f"No reference FASTA for condition(s): {missing}. "
            f"Provide --ref as a global default or --condition-ref NAME FASTA for each."
        )
    return condition_refs


def main():
    parser = build_parser()
    args = parser.parse_args()

    _DEFAULT_TITLE = "tRNA Alignment Pileup Heatmap"

    # Apply matplotlib styles before any plot is created
    _default_style = os.path.join(_PKG_DIR, 'plotting_styles', 'default.mplstyle')
    plt.style.use(_default_style)
    if getattr(args, 'mpl_style', None):
        plt.style.use(args.mpl_style)

    if args.command == 'run':
        conditions = pipeline.build_conditions(args.condition, parser.error)

        # Detect which conditions are pre-computed TSVs vs. BAM inputs.
        tsv_conditions = {n: paths[0] for n, paths in conditions.items()
                          if len(paths) == 1 and paths[0].endswith('.tsv')}
        bam_conditions = {n: paths for n, paths in conditions.items()
                          if n not in tsv_conditions}

        has_bam = bool(bam_conditions)
        has_tsv = bool(tsv_conditions)

        # Build per-condition ref map; validates --ref / --condition-ref coverage.
        condition_refs = _build_condition_refs(args, bam_conditions, parser)

        if has_bam and not (args.organism or args.cm or args.sprinzl_map):
            parser.error("--organism, --cm, or --sprinzl-map is required when any "
                         "--condition specifies BAM files.")

        # ── Collect Sprinzl-keyed rate dicts from all conditions ──────────────
        sprinzl_rates_by_condition = {}
        no_base_sets_by_condition  = {}
        counts_for_tsv             = None   # kept only for single BAM total-mode save
        stds_for_tsv               = None   # kept only for single BAM equal-mode save
        ref_to_sprinzl_for_tsv     = None   # ref_to_sprinzl for the condition that owns counts_for_tsv

        _sprinzl_cache       = {}  # {raw_ref_path: (axis, ref_to_sprinzl, mod_map)}
        per_cond_axes        = {}
        per_cond_r2s         = {}
        per_cond_mod         = {}
        per_cond_trimmed_seqs = {}
        _mapping_saved       = False

        if has_bam:
            print(f"Running pileup for {len(bam_conditions)} BAM condition(s) "
                  f"(merge-mode: {args.merge_mode})...")

        metric = 'match' if args.reference_match else 'mismatch'
        for cond_name, bam_paths in bam_conditions.items():
            cond_ref = condition_refs[cond_name]
            print(f"  Condition '{cond_name}': {len(bam_paths)} BAM(s)...")
            with _trim_context_for_ref(cond_ref, args, parser) as (ref_for_cmalign, trim5, trim3):
                if cond_ref not in _sprinzl_cache:
                    _sprinzl_cache[cond_ref] = _resolve_sprinzl_mapping(
                        args, ref_for_cmalign, parser)
                cond_axis, cond_r2s, cond_mod = _sprinzl_cache[cond_ref]
                per_cond_axes[cond_name]         = cond_axis
                per_cond_r2s[cond_name]          = cond_r2s
                per_cond_mod[cond_name]          = cond_mod
                per_cond_trimmed_seqs[cond_name] = pipeline.read_fasta_dict(ref_for_cmalign)

                if args.save_mapping and not _mapping_saved:
                    mapping_path = _resolve_output(args.save_mapping + '.tsv', args.outdir)
                    save_sprinzl_mapping(cond_axis, cond_r2s, mapping_path)
                    print(f"Saved Sprinzl mapping -> {mapping_path}")
                    _mapping_saved = True

                if args.merge_mode == 'equal':
                    r, s = pipeline.run_condition(
                        bam_paths, cond_ref, args.threads, 'equal',
                        min_q=args.min_q, metric=metric)
                    r = _filter_refs(pipeline.trim_arrays(r, trim5, trim3), args)
                    s = _filter_refs(pipeline.trim_arrays(s, trim5, trim3), args)
                    sprinzl_rates_by_condition[cond_name] = \
                        pipeline.project_to_sprinzl(r, cond_r2s)
                    if len(bam_conditions) == 1 and not has_tsv:
                        stds_for_tsv = pipeline.project_to_sprinzl(s, cond_r2s)
                else:
                    c = pipeline.run_condition(
                        bam_paths, cond_ref, args.threads, 'total',
                        min_q=args.min_q)
                    c = _filter_refs(pipeline.trim_arrays(c, trim5, trim3), args)
                    rate_fn = (pipeline.counts_to_accuracy if args.reference_match
                               else pipeline.counts_to_rates)
                    sprinzl_rates_by_condition[cond_name] = \
                        pipeline.project_to_sprinzl(rate_fn(c), cond_r2s)
                    if len(bam_conditions) == 1 and not has_tsv:
                        counts_for_tsv = c
                        ref_to_sprinzl_for_tsv = cond_r2s

        # ── TSV conditions: load directly ──────────────────────────────────
        tsv_axes = []
        for cond_name, tsv_path in tsv_conditions.items():
            print(f"  Condition '{cond_name}': loading from {tsv_path}...")
            sr, ax, nb, ref_order = heatmap.load_tsv(tsv_path)
            sr = _filter_refs(sr, args)
            nb = {k: v for k, v in nb.items() if k in sr}
            ref_order = [r for r in ref_order if r in sr]
            sprinzl_rates_by_condition[cond_name] = sr
            no_base_sets_by_condition[cond_name]  = nb
            tsv_axes.append(ax)

        # ── Build global Sprinzl axis (union across all sources) ───────────
        all_labels = set()
        for ax in per_cond_axes.values():
            all_labels.update(ax)
        for ax in tsv_axes:
            all_labels.update(ax)
        sprinzl_axis = sorted(all_labels, key=_sprinzl_sort_key)

        # ── no_base_sets for BAM conditions (recomputed against union axis) ─
        axis_set = set(sprinzl_axis)
        for cond_name in bam_conditions:
            no_base_sets_by_condition[cond_name] = {
                name: axis_set - set(labels)
                for name, labels in per_cond_r2s[cond_name].items()
            }

        # ── Auto-match refs by post-trim sequence identity ─────────────────
        per_cond_renames = pipeline.build_ref_name_map(per_cond_trimmed_seqs)
        for cond_name, renames in per_cond_renames.items():
            for orig, canonical in renames.items():
                if orig != canonical:
                    print(f"  Sequence match: '{cond_name}' '{orig}' → '{canonical}'")
        for cond_name in bam_conditions:
            rename = per_cond_renames[cond_name]
            sprinzl_rates_by_condition[cond_name] = {
                rename.get(k, k): v for k, v in sprinzl_rates_by_condition[cond_name].items()
            }
            no_base_sets_by_condition[cond_name] = {
                rename.get(k, k): v for k, v in no_base_sets_by_condition[cond_name].items()
            }
            if counts_for_tsv is not None:
                counts_for_tsv = {rename.get(k, k): v for k, v in counts_for_tsv.items()}
            if ref_to_sprinzl_for_tsv is not None:
                ref_to_sprinzl_for_tsv = {
                    rename.get(k, k): v for k, v in ref_to_sprinzl_for_tsv.items()
                }

        # ── mod_map: merge all per-condition maps ──────────────────────────
        mod_map = {}
        for cond_mod in per_cond_mod.values():
            for ref_name, label_map in cond_mod.items():
                mod_map.setdefault(ref_name, {}).update(label_map)

        # ── TSV save ───────────────────────────────────────────────────────
        if args.save_df:
            multi = len(conditions) > 1
            for cond_name in conditions:
                suffix = f"_{cond_name}" if multi else ""
                save_path = _resolve_output(args.save_df + suffix, args.outdir)
                ref_names_for_save = list(sprinzl_rates_by_condition[cond_name])
                no_base_for_save   = no_base_sets_by_condition.get(cond_name, {})
                if cond_name not in tsv_conditions and counts_for_tsv is not None:
                    heatmap.save_pileup(
                        counts_for_tsv, sprinzl_axis, ref_to_sprinzl_for_tsv, save_path)
                else:
                    heatmap.save_rates(
                        sprinzl_rates_by_condition[cond_name],
                        sprinzl_axis, ref_names_for_save, no_base_for_save,
                        save_path,
                        std_dict=stds_for_tsv if cond_name not in tsv_conditions else None)

            if multi:
                delta_prefix = _resolve_output(args.save_df, args.outdir)
                heatmap.save_deltas(
                    sprinzl_rates_by_condition, sprinzl_axis, delta_prefix,
                    no_base_sets_by_condition)

        output_path = _resolve_output(args.output, args.outdir)

        # ── Build unified ref_names (preserving order) ─────────────────────
        # Use order from first condition; later conditions may have extra refs
        # that get NaN in the matrix (same behaviour as mismatched BAMs).
        first_cond_rates = next(iter(sprinzl_rates_by_condition.values()))
        all_ref_names = list(first_cond_rates.keys())

        _plot_kwargs = dict(
            palette=args.palette, ylabel=args.ylabel,
            show_insertions=args.include_insertions,
            dpi=args.dpi, cell_size=args.cell_size, mod_map=mod_map,
            metric='match' if args.reference_match else 'mismatch',
            annotate=args.annotate,
            annotate_fontsize=args.annotate_fontsize,
            annotate_threshold=args.annotate_threshold,
        )

        if len(sprinzl_rates_by_condition) == 1:
            cond_name, sr = next(iter(sprinzl_rates_by_condition.items()))
            plot_title = cond_name if args.title == _DEFAULT_TITLE else args.title
            print(f"Generating heatmap -> {output_path}")
            heatmap.plot(sr, sprinzl_axis, all_ref_names,
                         no_base_sets_by_condition.get(cond_name, {}),
                         output_path, title=plot_title, **_plot_kwargs)
        else:
            if args.individual:
                base, ext = os.path.splitext(output_path)
                if not ext:
                    ext = '.pdf'
                for cond_name, sr in sprinzl_rates_by_condition.items():
                    ind_path = f"{base}_{cond_name}{ext}"
                    print(f"Generating individual heatmap -> {ind_path}")
                    ref_names_ind = list(sr.keys())
                    heatmap.plot(sr, sprinzl_axis, ref_names_ind,
                                 no_base_sets_by_condition.get(cond_name, {}),
                                 ind_path, title=cond_name, **_plot_kwargs)

            delta_title = args.title if args.title != _DEFAULT_TITLE else None
            print(f"Computing pairwise delta heatmaps "
                  f"for {len(sprinzl_rates_by_condition)} condition(s)...")
            heatmap.delta(
                sprinzl_rates_by_condition=sprinzl_rates_by_condition,
                sprinzl_axis=sprinzl_axis,
                output_prefix=output_path,
                no_base_sets_by_condition=no_base_sets_by_condition,
                title=delta_title,
                **_plot_kwargs,
            )

    elif args.command == 'inspect':
        with _adapter_trim_context(args, parser) as (ref_for_cmalign, _, _):
            sprinzl_axis, ref_to_sprinzl, mod_map = _resolve_sprinzl_mapping(
                args, ref_for_cmalign, parser)

            ref_to_sprinzl = _filter_refs(ref_to_sprinzl, args)
            sprinzl_axis = build_axis_from_mapping(ref_to_sprinzl)

            if args.save_mapping:
                mapping_path = _resolve_output(args.save_mapping + '.tsv', args.outdir)
                save_sprinzl_mapping(sprinzl_axis, ref_to_sprinzl, mapping_path)
                print(f"Saved Sprinzl mapping -> {mapping_path}")

            output_path = _resolve_output(args.output, args.outdir)
            title = args.title if args.title != _DEFAULT_TITLE else "tRNA Sprinzl Coverage"
            ref_to_seq = pipeline.read_fasta_dict(ref_for_cmalign) if args.include_sequence else None
            print(f"Generating Sprinzl coverage map -> {output_path}")
            heatmap.plot_sprinzl_coverage(
                sprinzl_axis, ref_to_sprinzl, output_path,
                palette=args.palette,
                ylabel=args.ylabel,
                title=title,
                show_insertions=args.include_insertions,
                ref_to_seq=ref_to_seq,
                seq_fontsize=args.seq_fontsize,
                dpi=args.dpi,
                cell_size=args.cell_size,
                mod_map=mod_map,
            )

    print("Done.")


if __name__ == "__main__":
    main()
