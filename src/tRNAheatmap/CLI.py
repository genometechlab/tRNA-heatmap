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
                                        load_sprinzl_mapping, build_axis_from_mapping)
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
    """Add --condition and --merge-mode flags to a subparser."""
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
        required=True,
        help="Path to reference tRNA FASTA file."
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
        help="Base name (no extension) to save pileup data as a .tsv file. "
             "Only valid with exactly one condition."
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


def _adapter_trim_context(args, parser):
    """
    CLI-side wrapper around pipeline.adapter_trimmed_ref.

    Validates adapter-related CLI flags using parser.error (so argparse
    formatting drives error messages), then returns the context manager that
    yields (effective_ref_path, trim_5, trim_3) and cleans up any temp file.
    """
    detect = getattr(args, 'detect_adapters', False)
    trim5  = getattr(args, 'trim_5', 0)
    trim3  = getattr(args, 'trim_3', 0)
    if detect and (trim5 or trim3):
        parser.error("--detect-adapters cannot be combined with --trim-5 or --trim-3.")
    if detect:
        seqs = pipeline.read_fasta(args.ref)
        if len(seqs) < 2:
            parser.error("--detect-adapters requires ≥2 sequences in the reference FASTA.")
    return pipeline.adapter_trimmed_ref(
        args.ref, trim_5=trim5, trim_3=trim3, detect=detect)


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
        if args.save_df and len(conditions) != 1:
            parser.error("--save-df is only valid with exactly one condition.")

        with _adapter_trim_context(args, parser) as (ref_for_cmalign, trim5, trim3):
            # Build {cond_name: rates_dict} uniformly; keep counts only when
            # we'll need them for the count-based TSV save (single-condition
            # total mode).
            rates_by_condition = {}
            counts_for_tsv     = None
            stds_for_tsv       = None
            print(f"Running pileup for {len(conditions)} condition(s) "
                  f"(merge-mode: {args.merge_mode})...")
            for cond_name, bam_paths in conditions.items():
                print(f"  Condition '{cond_name}': {len(bam_paths)} BAM(s)...")
                if args.merge_mode == 'equal':
                    r, s = pipeline.run_condition(bam_paths, args.ref, args.threads, 'equal')
                    rates_by_condition[cond_name] = _filter_refs(
                        pipeline.trim_arrays(r, trim5, trim3), args)
                    if len(conditions) == 1:
                        stds_for_tsv = _filter_refs(
                            pipeline.trim_arrays(s, trim5, trim3), args)
                else:
                    c = pipeline.run_condition(bam_paths, args.ref, args.threads, 'total')
                    c = _filter_refs(pipeline.trim_arrays(c, trim5, trim3), args)
                    rates_by_condition[cond_name] = pipeline.counts_to_rates(c)
                    if len(conditions) == 1:
                        counts_for_tsv = c

            sprinzl_axis, ref_to_sprinzl, mod_map = _resolve_sprinzl_mapping(
                args, ref_for_cmalign, parser)

            if args.save_mapping:
                mapping_path = _resolve_output(args.save_mapping + '.tsv', args.outdir)
                save_sprinzl_mapping(sprinzl_axis, ref_to_sprinzl, mapping_path)
                print(f"Saved Sprinzl mapping -> {mapping_path}")

            if args.save_df:
                save_path = _resolve_output(args.save_df, args.outdir)
                cond_name = next(iter(rates_by_condition))
                if counts_for_tsv is not None:
                    heatmap.save_pileup(counts_for_tsv, sprinzl_axis, ref_to_sprinzl, save_path)
                else:
                    heatmap.save_rates(rates_by_condition[cond_name], sprinzl_axis,
                                       ref_to_sprinzl, save_path, std_dict=stds_for_tsv)

            output_path = _resolve_output(args.output, args.outdir)

            if len(rates_by_condition) == 1:
                cond_name, rates = next(iter(rates_by_condition.items()))
                plot_title = cond_name if args.title == _DEFAULT_TITLE else args.title
                print(f"Generating heatmap -> {output_path}")
                heatmap.plot(rates, sprinzl_axis, ref_to_sprinzl, output_path,
                             palette=args.palette,
                             ylabel=args.ylabel,
                             title=plot_title,
                             show_insertions=args.include_insertions,
                             dpi=args.dpi, cell_size=args.cell_size, mod_map=mod_map)
            else:
                if args.individual:
                    base, ext = os.path.splitext(output_path)
                    if not ext:
                        ext = '.pdf'
                    for cond_name, rates in rates_by_condition.items():
                        ind_path = f"{base}_{cond_name}{ext}"
                        print(f"Generating individual heatmap -> {ind_path}")
                        heatmap.plot(rates, sprinzl_axis, ref_to_sprinzl, ind_path,
                                     palette=args.palette, ylabel=args.ylabel,
                                     title=cond_name, show_insertions=args.include_insertions,
                                     dpi=args.dpi, cell_size=args.cell_size, mod_map=mod_map)

                delta_title = args.title if args.title != _DEFAULT_TITLE else None
                print(f"Computing pairwise delta heatmaps for {len(rates_by_condition)} condition(s)...")
                heatmap.delta(
                    rates_by_condition=rates_by_condition,
                    sprinzl_axis=sprinzl_axis,
                    ref_to_sprinzl=ref_to_sprinzl,
                    output_prefix=output_path,
                    palette=args.palette,
                    ylabel=args.ylabel,
                    title=delta_title,
                    show_insertions=args.include_insertions,
                    dpi=args.dpi,
                    cell_size=args.cell_size,
                    mod_map=mod_map,
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
