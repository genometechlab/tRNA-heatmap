#!/usr/bin/env python3
"""
CLI entry point for tRNA pileup heatmap generation.

Subcommands
-----------
  run        Full pipeline: pileup → Sprinzl mapping → heatmap (+ optional save)
  plot       Plot heatmap from a pre-computed NPZ (skips pileup and cmalign)
  delta      Pairwise delta heatmaps from multiple pre-computed NPZ files
  run-delta  Full pipeline for ≥2 BAMs → pairwise delta heatmap in one step
  inspect    Dry-run: visualize Sprinzl coordinate coverage without any pileup

Examples
--------
  # Full pipeline with a bundled organism model
  python CLI.py run \\
      --ref ref.fa --bam aligned.bam --organism eukaryotic \\
      --output heatmap.png --save-df pileup

  # Delta heatmap directly from two BAMs (no intermediate NPZ files)
  python CLI.py run-delta \\
      --ref ref.fa --bam condA.bam condB.bam --organism eukaryotic \\
      --output results/delta

  # Plot only (no pileup or Infernal required)
  python CLI.py plot --load pileup.npz --output heatmap.png \\
      --palette dark-high --title "My Experiment" --ylabel "Isotype"

  # Delta heatmaps between all pairs of NPZ files
  python CLI.py delta --load A.npz B.npz C.npz --output results/delta

  # Inspect Sprinzl coordinate coverage (no BAM or pileup needed)
  python CLI.py inspect --ref ref.fa --organism eukaryotic --output coverage.png
"""

import argparse
import os
import tempfile
from .pileup_engine import pileup
from .calculate_tRNA_positions import (get_sprinzl_mapping, save_sprinzl_mapping,
                                        load_sprinzl_mapping, build_axis_from_mapping)
from . import heatmap

_PKG_DIR = os.path.dirname(os.path.abspath(__file__))
BUNDLED_MODELS = {
    'eukaryotic': os.path.join(_PKG_DIR, 'euk-num.cm'),
    'prokaryotic': os.path.join(_PKG_DIR, 'bact-num.cm'),
    'archaeal':    os.path.join(_PKG_DIR, 'arch-num.cm'),
}

PALETTE_CHOICES = list(heatmap.PALETTES.keys())


# ---------------------------------------------------------------------------
# FASTA / adapter utilities
# ---------------------------------------------------------------------------

def _read_fasta(ref_fasta):
    """Parse FASTA file. Returns list of (name, seq) tuples."""
    seqs, name, buf = [], None, []
    with open(ref_fasta) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if name:
                    seqs.append((name, ''.join(buf)))
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line)
    if name:
        seqs.append((name, ''.join(buf)))
    return seqs


def _detect_adapters(sequences):
    """
    Find common 5' prefix and 3' suffix lengths across all sequences.

    The terminal 3' CCA is always preserved: if the detected suffix ends in
    CCA, trim3 is reduced by 3 so the biologically required CCA tail is kept.
    """
    seqs = [s for _, s in sequences]
    trim5 = 0
    for bases in zip(*seqs):
        if len(set(bases)) == 1:
            trim5 += 1
        else:
            break
    trim3 = 0
    for bases in zip(*[s[::-1] for s in seqs]):
        if len(set(bases)) == 1:
            trim3 += 1
        else:
            break
    # Protect the CCA tail: if the detected suffix begins with CCA (reading
    # 5'→3'), it is the biologically required mature tRNA tail, not adapter.
    # e.g. suffix = 'CCAGGCTTC' → trim only 'GGCTTC'; suffix = 'CCA' → trim nothing.
    if trim3 >= 3:
        L = len(seqs[0])
        if seqs[0][L - trim3 : L - trim3 + 3].upper() == 'CCA':
            trim3 -= 3
    return trim5, trim3


def _write_trimmed_fasta(ref_fasta, trim5, trim3):
    """Write trim-adjusted sequences to a temp file. Caller must os.unlink()."""
    seqs = _read_fasta(ref_fasta)
    fd, tmp_path = tempfile.mkstemp(suffix='.fa')
    with os.fdopen(fd, 'w') as f:
        for name, seq in seqs:
            end = len(seq) - trim3 if trim3 > 0 else len(seq)
            f.write(f'>{name}\n{seq[trim5:end]}\n')
    return tmp_path


def _trim_arrays(arrays, trim5, trim3):
    """Slice pileup arrays to remove adapter positions."""
    if trim5 == 0 and trim3 == 0:
        return arrays
    result = {}
    for name, arr in arrays.items():
        end = arr.shape[1] - trim3 if trim3 > 0 else arr.shape[1]
        result[name] = arr[:, trim5:end]
    return result


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
             "They are always stored in NPZ/TSV files."
    )
    p.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Output image resolution in DPI. Default: 300."
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

    if has_map:
        tsv_mapping = load_sprinzl_mapping(args.sprinzl_map)
        if not has_cm:
            ref_to_sprinzl = tsv_mapping
            sprinzl_axis = build_axis_from_mapping(ref_to_sprinzl)
        else:
            # Patch mode: TSV entries overwrite cmalign results for listed refs
            ref_to_sprinzl.update(tsv_mapping)

    return sprinzl_axis, ref_to_sprinzl


def build_parser():
    parser = argparse.ArgumentParser(
        description="tRNA pileup heatmap tool.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest='command', required=True)

    # ------------------------------------------------------------------
    # 'run' subcommand — full pipeline
    # ------------------------------------------------------------------
    run_p = subparsers.add_parser(
        'run',
        help='Full pipeline: pileup + Sprinzl mapping + heatmap.',
    )
    run_p.add_argument(
        "--ref", "-r",
        required=True,
        help="Path to reference tRNA FASTA file."
    )
    run_p.add_argument(
        "--bam", "-b",
        required=True,
        help="Path to aligned BAM file (must be sorted)."
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
        default="heatmap.png",
        help="Output heatmap file. Format from extension (.png, .pdf, .svg). "
             "Default: heatmap.png"
    )
    run_p.add_argument(
        "--save-df", "-s",
        default=None,
        dest="save_df",
        help="Base name (no extension) to save pileup as both .npz and .tsv."
    )
    run_p.add_argument(
        "--threads", "-t",
        type=int,
        default=1,
        help="Threads for pileup engine. Default: 1"
    )
    _add_adapter_flags(run_p)
    _add_sprinzl_override_flags(run_p)
    _add_plot_flags(run_p)

    # ------------------------------------------------------------------
    # 'plot' subcommand — plot from pre-computed NPZ
    # ------------------------------------------------------------------
    plot_p = subparsers.add_parser(
        'plot',
        help='Plot heatmap from a pre-computed NPZ file.',
    )
    plot_p.add_argument(
        "--load", "-l",
        required=True,
        help="Path to a pre-computed .npz file from 'run --save-df'."
    )
    plot_p.add_argument(
        "--output", "-o",
        default="heatmap.png",
        help="Output heatmap file. Default: heatmap.png"
    )
    _add_plot_flags(plot_p)

    # ------------------------------------------------------------------
    # 'delta' subcommand — pairwise delta heatmaps
    # ------------------------------------------------------------------
    delta_p = subparsers.add_parser(
        'delta',
        help='Pairwise delta heatmaps across multiple pre-computed NPZ files.',
    )
    delta_p.add_argument(
        "--load", "-l",
        nargs='+',
        required=True,
        metavar="NPZ",
        help="Two or more pre-computed .npz files."
    )
    delta_p.add_argument(
        "--output", "-o",
        default="delta",
        help="Output prefix (with optional extension). "
             "Produces {prefix}_{A}_vs_{B}.{ext} per pair. Default: 'delta'."
    )
    _add_plot_flags(delta_p)

    # ------------------------------------------------------------------
    # 'run-delta' subcommand — full pipeline straight to delta heatmap
    # ------------------------------------------------------------------
    rundelta_p = subparsers.add_parser(
        'run-delta',
        help='Full pipeline for ≥2 BAMs, producing pairwise delta heatmaps directly.',
    )
    rundelta_p.add_argument(
        "--ref", "-r",
        required=True,
        help="Path to reference tRNA FASTA file (shared across all BAMs)."
    )
    rundelta_p.add_argument(
        "--bam", "-b",
        nargs='+',
        required=True,
        metavar="BAM",
        help="Two or more sorted BAM files to compare."
    )
    cm_group_rd = rundelta_p.add_mutually_exclusive_group(required=False)
    cm_group_rd.add_argument(
        "--organism", "-g",
        choices=["eukaryotic", "prokaryotic", "archaeal"],
        dest="organism",
        help="Use a bundled covariance model: eukaryotic (euk-num.cm), "
             "prokaryotic (bact-num.cm), or archaeal (arch-num.cm)."
    )
    cm_group_rd.add_argument(
        "--cm", "-c",
        dest="cm",
        metavar="FILE",
        help="Path to a custom Infernal covariance model (.cm file)."
    )
    rundelta_p.add_argument(
        "--output", "-o",
        default="delta",
        help="Output prefix (with optional extension). "
             "Produces {prefix}_{stemA}_vs_{stemB}.{ext} per pair. Default: 'delta'."
    )
    rundelta_p.add_argument(
        "--threads", "-t",
        type=int,
        default=1,
        help="Threads for pileup engine. Default: 1"
    )
    _add_adapter_flags(rundelta_p)
    _add_sprinzl_override_flags(rundelta_p)
    _add_plot_flags(rundelta_p)

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
        default="sprinzl_coverage.png",
        help="Output heatmap file. Format from extension (.png, .pdf, .svg). "
             "Default: sprinzl_coverage.png"
    )
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


def _apply_adapter_trim(args, parser):
    """
    Resolve adapter trimming from CLI flags.

    Returns (trim5, trim3, tmp_ref_path_or_None). If a temp file is returned,
    the caller must os.unlink() it in a finally block.
    """
    detect = getattr(args, 'detect_adapters', False)
    trim5  = getattr(args, 'trim_5', 0)
    trim3  = getattr(args, 'trim_3', 0)

    if detect and (trim5 or trim3):
        parser.error("--detect-adapters cannot be combined with --trim-5 or --trim-3.")

    if detect:
        seqs = _read_fasta(args.ref)
        if len(seqs) < 2:
            parser.error("--detect-adapters requires ≥2 sequences in the reference FASTA.")
        trim5, trim3 = _detect_adapters(seqs)
        print(f"Detected adapters: 5' trim = {trim5} bp, 3' trim = {trim3} bp")

    if trim5 == 0 and trim3 == 0:
        return 0, 0, None

    tmp = _write_trimmed_fasta(args.ref, trim5, trim3)
    return trim5, trim3, tmp


def main():
    parser = build_parser()
    args = parser.parse_args()

    if args.command == 'run':
        trim5, trim3, tmp_ref = _apply_adapter_trim(args, parser)
        ref_for_cmalign = tmp_ref or args.ref
        try:
            print(f"Running pileup on {args.bam} against {args.ref} "
                  f"using {args.threads} thread(s)...")
            arrays = pileup(args.bam, args.ref, args.threads)
            arrays = _trim_arrays(arrays, trim5, trim3)

            sprinzl_axis, ref_to_sprinzl = _resolve_sprinzl_mapping(
                args, ref_for_cmalign, parser)

            if args.save_mapping:
                mapping_path = _resolve_output(args.save_mapping + '.tsv', args.outdir)
                save_sprinzl_mapping(sprinzl_axis, ref_to_sprinzl, mapping_path)
                print(f"Saved Sprinzl mapping -> {mapping_path}")

            if args.save_df:
                save_path = _resolve_output(args.save_df, args.outdir)
                print(f"Saving pileup data -> {save_path}.npz + {save_path}.tsv")
                heatmap.save_pileup(arrays, sprinzl_axis, ref_to_sprinzl, save_path)

            output_path = _resolve_output(args.output, args.outdir)
            print(f"Generating heatmap -> {output_path}")
            heatmap.plot(arrays, sprinzl_axis, ref_to_sprinzl, output_path,
                         palette=args.palette,
                         ylabel=args.ylabel,
                         title=args.title,
                         show_insertions=args.include_insertions,
                         dpi=args.dpi)
        finally:
            if tmp_ref:
                os.unlink(tmp_ref)

    elif args.command == 'plot':
        output_path = _resolve_output(args.output, args.outdir)
        print(f"Loading pre-computed pileup from {args.load}...")
        heatmap.plot_from_npz(args.load, output_path,
                              palette=args.palette,
                              ylabel=args.ylabel,
                              title=args.title,
                              show_insertions=args.include_insertions,
                              dpi=args.dpi)

    elif args.command == 'delta':
        if len(args.load) < 2:
            parser.error("delta requires at least 2 NPZ files via --load.")
        output_prefix = _resolve_output(args.output, args.outdir)
        print(f"Computing pairwise delta heatmaps for {len(args.load)} file(s)...")
        heatmap.delta_heatmap(
            npz_paths=args.load,
            output_prefix=output_prefix,
            palette=args.palette,
            ylabel=args.ylabel,
            title=args.title if args.title != "tRNA Alignment Pileup Heatmap" else None,
            show_insertions=args.include_insertions,
            dpi=args.dpi,
        )

    elif args.command == 'run-delta':
        if len(args.bam) < 2:
            parser.error("run-delta requires at least 2 BAM files via --bam.")
        trim5, trim3, tmp_ref = _apply_adapter_trim(args, parser)
        ref_for_cmalign = tmp_ref or args.ref
        try:
            sprinzl_axis, ref_to_sprinzl = _resolve_sprinzl_mapping(
                args, ref_for_cmalign, parser)

            if args.save_mapping:
                mapping_path = _resolve_output(args.save_mapping + '.tsv', args.outdir)
                save_sprinzl_mapping(sprinzl_axis, ref_to_sprinzl, mapping_path)
                print(f"Saved Sprinzl mapping -> {mapping_path}")

            arrays_by_stem = {}
            for bam_path in args.bam:
                stem = os.path.splitext(os.path.basename(bam_path))[0]
                print(f"Running pileup on {bam_path} using {args.threads} thread(s)...")
                arr = pileup(bam_path, args.ref, args.threads)
                arrays_by_stem[stem] = _trim_arrays(arr, trim5, trim3)

            output_prefix = _resolve_output(args.output, args.outdir)
            print(f"Computing pairwise delta heatmaps for {len(args.bam)} BAM(s)...")
            heatmap.delta_from_arrays(
                arrays_by_stem=arrays_by_stem,
                sprinzl_axis=sprinzl_axis,
                ref_to_sprinzl=ref_to_sprinzl,
                output_prefix=output_prefix,
                palette=args.palette,
                ylabel=args.ylabel,
                title=args.title if args.title != "tRNA Alignment Pileup Heatmap" else None,
                show_insertions=args.include_insertions,
                dpi=args.dpi,
            )
        finally:
            if tmp_ref:
                os.unlink(tmp_ref)

    elif args.command == 'inspect':
        trim5, trim3, tmp_ref = _apply_adapter_trim(args, parser)
        ref_for_cmalign = tmp_ref or args.ref
        try:
            sprinzl_axis, ref_to_sprinzl = _resolve_sprinzl_mapping(
                args, ref_for_cmalign, parser)

            if args.save_mapping:
                mapping_path = _resolve_output(args.save_mapping + '.tsv', args.outdir)
                save_sprinzl_mapping(sprinzl_axis, ref_to_sprinzl, mapping_path)
                print(f"Saved Sprinzl mapping -> {mapping_path}")

            output_path = _resolve_output(args.output, args.outdir)
            title = args.title if args.title != "tRNA Alignment Pileup Heatmap" else "tRNA Sprinzl Coverage"
            print(f"Generating Sprinzl coverage map -> {output_path}")
            heatmap.plot_sprinzl_coverage(
                sprinzl_axis, ref_to_sprinzl, output_path,
                palette=args.palette,
                ylabel=args.ylabel,
                title=title,
                show_insertions=args.include_insertions,
                dpi=args.dpi,
            )
        finally:
            if tmp_ref:
                os.unlink(tmp_ref)

    print("Done.")


if __name__ == "__main__":
    main()
