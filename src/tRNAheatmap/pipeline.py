"""
Pure pipeline data functions.

Anything that touches reads, references, conditions, or merging — but does NOT
parse argv, build matplotlib figures, or know about argparse — lives here.
CLI.py imports from this module; tests can import these helpers directly
without going through subprocess.

Public API
----------
  read_fasta(path)              → list[(name, seq)]
  read_fasta_dict(path)         → dict[name, seq]
  detect_adapters(seqs)         → (trim_5, trim_3) protecting the CCA tail
  adapter_trimmed_ref(...)      → context manager yielding (ref_path, t5, t3)
  trim_arrays(arrays, t5, t3)   → slice 1D rate arrays or 2D count arrays
  apply_ref_filter(arrays, ...) → keep/drop refs by include/exclude file
  merge_total(arrays_list)      → element-wise sum of count arrays
  merge_equal(arrays_list)      → mean of per-BAM rates with std dev
  counts_to_rates(arrays)       → derive 1D mismatch rate from 4×L count arrays
  run_condition(bams, ref, threads, merge_mode) → pileup each BAM + merge
  build_conditions(condition_groups, error) → parse repeated --condition groups
"""

import contextlib
import os
import tempfile

import numpy as np

from .pileup_engine import pileup


# ---------------------------------------------------------------------------
# FASTA / adapter utilities
# ---------------------------------------------------------------------------

def read_fasta(path):
    """Parse FASTA file into a list of (name, seq) tuples."""
    seqs, name, buf = [], None, []
    with open(path) as fh:
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


def read_fasta_dict(path):
    """Return {name: sequence} for every record in a FASTA file."""
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line)
    return {k: ''.join(v) for k, v in seqs.items()}


def detect_adapters(sequences):
    """
    Find common 5' prefix and 3' suffix lengths across all sequences.

    The terminal 3' CCA is always preserved: if the detected suffix ends in
    CCA, trim_3 is reduced by 3 so the biologically required CCA tail is kept.
    """
    seqs = [s for _, s in sequences]
    trim_5 = 0
    for bases in zip(*seqs):
        if len(set(bases)) == 1:
            trim_5 += 1
        else:
            break
    trim_3 = 0
    for bases in zip(*[s[::-1] for s in seqs]):
        if len(set(bases)) == 1:
            trim_3 += 1
        else:
            break
    # Protect the CCA tail: if the detected suffix begins with CCA (reading
    # 5'→3'), it is the biologically required mature tRNA tail, not adapter.
    if trim_3 >= 3:
        L = len(seqs[0])
        if seqs[0][L - trim_3 : L - trim_3 + 3].upper() == 'CCA':
            trim_3 -= 3
    return trim_5, trim_3


def _write_trimmed_fasta(ref_fasta, trim_5, trim_3):
    """Write trim-adjusted sequences to a temp file. Caller must os.unlink()."""
    seqs = read_fasta(ref_fasta)
    fd, tmp_path = tempfile.mkstemp(suffix='.fa')
    with os.fdopen(fd, 'w') as f:
        for name, seq in seqs:
            end = len(seq) - trim_3 if trim_3 > 0 else len(seq)
            f.write(f'>{name}\n{seq[trim_5:end]}\n')
    return tmp_path


@contextlib.contextmanager
def adapter_trimmed_ref(ref_fasta, *, trim_5=0, trim_3=0, detect=False):
    """
    Resolve adapter trimming and yield (effective_ref_path, trim_5, trim_3).

    If any trimming is applied the yielded path points at a temp file that is
    deleted when the context exits. If no trimming is needed the original path
    is yielded unchanged.

    Raises ValueError if `detect` is combined with explicit trim_5/trim_3 or
    if `detect` is requested but the FASTA has fewer than 2 sequences.
    """
    if detect and (trim_5 or trim_3):
        raise ValueError("detect_adapters cannot be combined with trim_5 or trim_3.")

    if detect:
        seqs = read_fasta(ref_fasta)
        if len(seqs) < 2:
            raise ValueError("detect_adapters requires ≥2 sequences in the reference FASTA.")
        trim_5, trim_3 = detect_adapters(seqs)
        print(f"Detected adapters: 5' trim = {trim_5} bp, 3' trim = {trim_3} bp")

    if trim_5 == 0 and trim_3 == 0:
        yield ref_fasta, 0, 0
        return

    tmp = _write_trimmed_fasta(ref_fasta, trim_5, trim_3)
    try:
        yield tmp, trim_5, trim_3
    finally:
        os.unlink(tmp)


def trim_arrays(arrays, trim_5, trim_3):
    """
    Slice pileup arrays to remove adapter positions.

    Handles both 2D count arrays (4, L) and 1D rate arrays (L,) produced
    by equal-weight merging.
    """
    if trim_5 == 0 and trim_3 == 0:
        return arrays
    out = {}
    for name, arr in arrays.items():
        if arr.ndim == 1:
            end = len(arr) - trim_3 if trim_3 > 0 else len(arr)
            out[name] = arr[trim_5:end]
        else:
            end = arr.shape[1] - trim_3 if trim_3 > 0 else arr.shape[1]
            out[name] = arr[:, trim_5:end]
    return out


# ---------------------------------------------------------------------------
# Reference filtering
# ---------------------------------------------------------------------------

def apply_ref_filter(arrays, include_file=None, exclude_file=None):
    """
    Keep or drop reference names from an arrays dict.

    Works for both 4-row count arrays and 1D rate arrays — only the dict keys
    are touched. Returns arrays unchanged if neither file is given.
    `include_file` preserves the file's listed order. `exclude_file` preserves
    the dict's order minus the listed names.
    """
    if not include_file and not exclude_file:
        return arrays
    if include_file:
        with open(include_file) as fh:
            ordered = [ln.strip() for ln in fh if ln.strip() and not ln.startswith('#')]
        return {k: arrays[k] for k in ordered if k in arrays}
    with open(exclude_file) as fh:
        drop = {ln.strip() for ln in fh if ln.strip() and not ln.startswith('#')}
    return {k: v for k, v in arrays.items() if k not in drop}


# ---------------------------------------------------------------------------
# Condition merging
# ---------------------------------------------------------------------------

def merge_total(arrays_list):
    """Sum raw pileup count arrays across a list of BAM results."""
    merged = {name: arr.copy() for name, arr in arrays_list[0].items()}
    for a in arrays_list[1:]:
        for name in merged:
            merged[name] += a[name]
    return merged


def merge_equal(arrays_list):
    """
    Average per-BAM mismatch rates per position.

    Returns (means, stds) — two dicts mapping ref_name → 1D float64 array.
    NaN where a position had no coverage in any replicate.
    stds is NaN (not 0) for single-BAM conditions (undefined std dev).
    """
    means, stds = {}, {}
    for name in arrays_list[0]:
        per_bam = []
        for a in arrays_list:
            arr   = a[name]
            match = arr[0].astype(np.float64)
            mm    = arr[1].astype(np.float64)
            total = match + mm
            with np.errstate(invalid='ignore', divide='ignore'):
                per_bam.append(np.where(total > 0, mm / total, np.nan))
        with np.errstate(all='ignore'):
            means[name] = np.nanmean(per_bam, axis=0)
            stds[name]  = (np.nanstd(per_bam, axis=0) if len(per_bam) > 1
                           else np.full_like(means[name], np.nan))
    return means, stds


def counts_to_rates(count_arrays):
    """
    Convert {name: (4, L) count array} to {name: (L,) mismatch-rate array}.

    Both 'total' and 'equal' merge modes converge on this rate representation
    before plotting; equal mode produces it directly, total mode derives it
    from merged raw counts. NaN where a position had no coverage.
    """
    rates = {}
    for name, arr in count_arrays.items():
        match = arr[0].astype(np.float64)
        mm    = arr[1].astype(np.float64)
        total = match + mm
        with np.errstate(invalid='ignore', divide='ignore'):
            rates[name] = np.where(total > 0, mm / total, np.nan)
    return rates


def run_condition(bam_paths, ref, threads, merge_mode):
    """Run pileup for each BAM in a condition, then merge."""
    arrays_list = [pileup(b, ref, threads) for b in bam_paths]
    if merge_mode == 'equal':
        return merge_equal(arrays_list)
    return merge_total(arrays_list)


# ---------------------------------------------------------------------------
# Condition input normalization
# ---------------------------------------------------------------------------

def build_conditions(condition_groups, error):
    """
    Parse repeated --condition groups into an ordered dict {cond_name: [bam_paths]}.

    Parameters
    ----------
    condition_groups : list[list[str]] or None
        argparse `--condition` action='append' result; each inner list is
        [name, bam_path, bam_path, ...].
    error : callable
        Called as error(msg) on validation failure (typically parser.error).

    Returns
    -------
    OrderedDict-style dict[str, list[str]]
    """
    if not condition_groups:
        error("Provide at least one --condition NAME BAM [BAM ...] group.")
    conds = {}
    for group in condition_groups:
        if len(group) < 2:
            error(
                f"--condition requires a condition name followed by at least "
                f"one BAM file. Got: {group!r}"
            )
        name, bams = group[0], group[1:]
        if name in conds:
            error(f"Duplicate condition name: {name!r}")
        conds[name] = bams
    return conds
