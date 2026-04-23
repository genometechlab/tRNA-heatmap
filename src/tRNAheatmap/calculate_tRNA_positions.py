#TODO: Add a flag to either include or exclude insertions from the infernal alignment.
"""
Sprinzl coordinate assignment for tRNA sequences.

This module uses components from tRNAviz (UCSC Lowe Lab) for mapping
tRNA nucleotides to canonical Sprinzl positional numbering via Infernal
covariance model alignment.

The following files are redistributed from tRNAviz under LGPL-3.0:
  - tRNA_position.py : Sprinzl position mapping from SS_cons annotations
  - euk-num.cm       : Eukaryotic tRNA numbering covariance model
  - bact-num.cm      : Bacterial tRNA numbering covariance model
  - arch-num.cm      : Archaeal tRNA numbering covariance model

Source: https://github.com/UCSC-LoweLab/tRNAviz
License: LGPL-3.0 (https://github.com/UCSC-LoweLab/tRNAviz/blob/master/LICENSE)

Citation:
  Lin BY, Chan PP, Lowe TM. (2019) tRNAviz: explore and visualize tRNA
  sequence features. Nucleic Acids Res. 47(W1):W542-W547.
  https://doi.org/10.1093/nar/gkz438

External dependency:
  Infernal (cmalign) - http://eddylab.org/infernal/
  Nawrocki EP, Eddy SR. (2013) Infernal 1.1: 100-fold faster RNA
  homology searches. Bioinformatics. 29:2933-2935.

Public API
----------
  get_sprinzl_mapping(ref_fasta, cm_model)
      Run cmalign, parse SS_cons, return the Sprinzl x-axis and a per-reference
      mapping from physical reference position to Sprinzl label.

      Returns
      -------
      sprinzl_axis : list[str]
          Ordered Sprinzl labels for every consensus column (the shared x-axis).
      ref_to_sprinzl : dict[str, list[str]]
          {ref_name: [sprinzl_label, ...]} — one label per physical reference
          position (length == pileup array width for that reference).

  save_sprinzl_mapping(sprinzl_axis, ref_to_sprinzl, path)
      Write the mapping to a human-readable TSV for inspection or reuse.

  load_sprinzl_mapping(path)
      Load a TSV written by save_sprinzl_mapping (or hand-edited by a user).
      Returns (sprinzl_axis_or_None, ref_to_sprinzl).
"""

import re
import subprocess
import tempfile
import os
from collections import defaultdict


# ---------------------------------------------------------------------------
# Classes (from tRNAviz)
# ---------------------------------------------------------------------------

class Position:
  def __init__(self, position, sprinzl, index=-1, paired=False):
    self.position = position
    self.sprinzl = sprinzl
    self.index = index
    self.paired = paired

  def __str__(self):
    return "Position {} (Sprinzl: {})".format(self.position, self.sprinzl)


class Region:
  def __init__(self, lower, upper, name):
    self.lower = lower
    self.upper = upper
    self.name = name

  def __str__(self):
    return "Region {} ({} ~ {})".format(self.name, self.lower, self.upper)


# ---------------------------------------------------------------------------
# annotate_positions — from tRNAviz
# ---------------------------------------------------------------------------

def annotate_positions(ss):
  loop_indices = []
  for i, indices in enumerate([r.span() for r in re.finditer(r'([\(\.]+\()|(<[<\.]+<)|[_\.]+|>[>\.]+>|\)[\)\.]+', ss)]):
    left, right = indices
    if re.search(r'[\(<_>\)]', ss[left:right]): loop_indices.append(indices)
  regions = ['acceptor', 'dstem', 'dloop', 'dstem', 'acstem', 'acloop', 'acstem', 'vstem', 'vloop', 'vstem', 'tpcstem', 'tpcloop', 'tpcstem', 'acceptor']
  regions = [Region(indices[0], indices[1], name) for indices, name in zip(loop_indices, regions[:])]
  region = regions[0]
  positions = []
  region_index = 0
  region_numbering = 0
  insert_index = 0
  sprinzl_positions = ['1:72', '2:71', '3:70', '4:69', '5:68', '6:67', '7:66', '8', '9', '10:25', '11:24', '12:23', '13:22', '14', '15', '16', '17', '17a', '18', '19', '20', '20a', '20b', '21', '22:13', '23:12', '24:11', '25:10', '26', '27:43', '28:42', '29:41', '30:40', '31:39', '32', '33', '34', '35', '36', '37', '38', '39:31', '40:30', '41:29', '42:28', '43:27', '44', '45', 'V11:V21', 'V12:V22', 'V13:V23', 'V14:V24', 'V15:V25', 'V16:V26', 'V17:V27', 'V1', 'V2', 'V3', 'V4', 'V5', 'V27:V17', 'V26:V16', 'V25:V15', 'V24:V14', 'V23:V13', 'V22:V12', 'V21:V11', '46', '47', '48', '49:65', '50:64', '51:63', '52:62', '53:61', '54', '55', '56', '57', '58', '59', '60', '61:53', '62:52', '63:51', '64:50', '65:49', '66:7', '67:6', '68:5', '69:4', '70:3', '71:2', '72:1', '73', '74', '75', '76']
  sprinzl_index = 0
  sprinzl_insert_index = 0
  reverse_bp_index = 1
  for position in range(len(ss)):
    if region_index < len(regions):
      region = regions[region_index]
    if ss[position] == '.':
      insert_index += 1
      sprinzl_insert_index += 1
      sprinzl = '{}i{}'.format(sprinzl_positions[sprinzl_index - 1].split(':')[0], sprinzl_insert_index)
      positions.append(Position(position=str(position + 1), sprinzl=sprinzl, index=len(positions), paired=False))
    else:
      if position < region.lower:
        positions.append(Position(position=str(position + 1), sprinzl=sprinzl_positions[sprinzl_index], index=len(positions), paired=False))
      elif position == region.lower:
        insert_index = 0
        if ss[position] == "(":
          while ss[regions[-1].upper - reverse_bp_index] == ".": reverse_bp_index += 1
          paired_base = regions[-1].upper - reverse_bp_index
          positions.append(Position(position='{}:{}'.format(position + 1, paired_base + 1), sprinzl=sprinzl_positions[sprinzl_index], index=len(positions), paired=True))
          region_numbering += 1
        elif ss[position] == "<":
          while ss[regions[region_index + 2].upper - reverse_bp_index] == '.': reverse_bp_index += 1
          paired_base = regions[region_index + 2].upper - reverse_bp_index
          positions.append(Position(position='{}:{}'.format(position + 1, paired_base + 1), sprinzl=sprinzl_positions[sprinzl_index], index=len(positions), paired=True))
          region_numbering += 1
        elif ss[position] in [')', '>']:
          pass
        else:
          positions.append(Position(position=str(position + 1), sprinzl=sprinzl_positions[sprinzl_index], index=len(positions), paired=False))
      elif position > region.lower and position <= region.upper - 1:
        if ss[position] == "(":
          while ss[regions[-1].upper - region_numbering - reverse_bp_index] == ".": reverse_bp_index += 1
          paired_base = regions[-1].upper - region_numbering - reverse_bp_index
          positions.append(Position(position='{}:{}'.format(position + 1, paired_base + 1), sprinzl=sprinzl_positions[sprinzl_index], index=len(positions), paired=True))
          region_numbering += 1
        elif ss[position] == "<":
          while ss[regions[region_index + 2].upper - region_numbering - reverse_bp_index] == ".": reverse_bp_index += 1
          paired_base = regions[region_index + 2].upper - region_numbering - reverse_bp_index
          positions.append(Position(position='{}:{}'.format(position + 1, paired_base + 1), sprinzl=sprinzl_positions[sprinzl_index], index=len(positions), paired=True))
          region_numbering += 1
        elif ss[position] in [')', '>']:
          pass
        else:
          positions.append(Position(position=str(position + 1), sprinzl=sprinzl_positions[sprinzl_index], index=len(positions), paired=False))
      else:
        positions.append(Position(position=position + 1, sprinzl=sprinzl_positions[sprinzl_index], index=len(positions), paired=False))
      sprinzl_index += 1
      sprinzl_insert_index = 0
    if position == region.upper - 1:
      region_index += 1
      reverse_bp_index = 1
      region_numbering = 0
      insert_index = 0
  return positions


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def parse_alignment_cons(alignment_file):
  """Read SS_cons lines from a Stockholm alignment and return (ss, positions)."""
  ss = ''
  with open(alignment_file) as fh:
    for line in fh:
      if line[0:12] == '#=GC SS_cons':
        ss += line.strip().split()[-1]
  return ss, annotate_positions(ss)


def build_col_to_sprinzl(ss, positions):
  """
  Build a list mapping each alignment column index to its Sprinzl label.

  The SS_cons string has one character per alignment column. The Position list
  produced by annotate_positions() skips ')' and '>' characters (closing sides
  of base pairs) — those columns are assigned the second half of their partner's
  Sprinzl label (e.g., column for ')' paired with '(' at Sprinzl '1:72' gets
  label '72').

  Parameters
  ----------
  ss : str
      SS_cons string (full length, one char per alignment column).
  positions : list[Position]
      Output of annotate_positions(ss).

  Returns
  -------
  list[str] of length len(ss)
  """
  col_to_sprinzl = [''] * len(ss)

  # First pass: assign labels for all non-skip columns (not ')' or '>')
  pos_iter = iter(positions)
  for col, char in enumerate(ss):
    if char not in (')', '>'):
      p = next(pos_iter)
      # For paired positions the sprinzl label is like '1:72'; use first part
      col_to_sprinzl[col] = p.sprinzl.split(':')[0]

  # Second pass: fill in ')' and '>' columns using paired Position partner info
  for p in positions:
    if p.paired:
      # p.position is '1:72' (1-indexed); partner column is the second number
      partner_col = int(p.position.split(':')[1]) - 1
      col_to_sprinzl[partner_col] = p.sprinzl.split(':')[1]

  return col_to_sprinzl


def build_sprinzl_mapping(alignment_file, col_to_sprinzl):
  """
  Read aligned sequences from a Stockholm file and build per-reference
  Sprinzl label lists.

  Parameters
  ----------
  alignment_file : str
      Path to the .sto file produced by cmalign.
  col_to_sprinzl : list[str]
      Output of build_col_to_sprinzl(); one label per alignment column.

  Returns
  -------
  sprinzl_axis : list[str]
      Ordered Sprinzl labels for every consensus column (empty-string entries
      excluded — these are positions that are gap-only in the consensus).
  ref_to_sprinzl : dict[str, list[str]]
      {seq_name: [sprinzl_label, ...]} where the list has one entry per
      physical reference position (i.e., same length as the pileup array
      for that reference).
  """
  seqs = defaultdict(str)
  seqnames = []

  with open(alignment_file) as fh:
    for line in fh:
      if line[0] in ('#', '\n', '/'):
        continue
      parts = line.strip().split()
      if len(parts) != 2:
        continue
      seqname, seq = parts
      if seqname not in seqnames:
        seqnames.append(seqname)
      seqs[seqname] += seq

  # Consensus x-axis: labels for columns that are not purely gap ('-' in all seqs)
  sprinzl_axis = [label for label in col_to_sprinzl if label != '']

  # Per-reference mapping: one Sprinzl label per non-gap character in that sequence
  ref_to_sprinzl = {}
  for seqname in seqnames:
    seq = seqs[seqname]
    labels = [
      col_to_sprinzl[col]
      for col, char in enumerate(seq)
      if char not in ('-', '.')
    ]
    ref_to_sprinzl[seqname] = labels

  return sprinzl_axis, ref_to_sprinzl


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def get_sprinzl_mapping(ref_fasta, cm_model):
  """
  Align ref_fasta to cm_model with cmalign, parse the consensus secondary
  structure, and return the Sprinzl x-axis and per-reference position mapping.

  Parameters
  ----------
  ref_fasta : str
      Path to the reference tRNA FASTA file.
  cm_model : str
      Path to an Infernal covariance model (.cm) optimised for tRNA numbering.

  Returns
  -------
  sprinzl_axis : list[str]
      Ordered Sprinzl labels for every consensus column (the shared x-axis).
  ref_to_sprinzl : dict[str, list[str]]
      {ref_name: [sprinzl_label, ...]} — one label per physical reference
      position (length == pileup array width for that reference).
  """
  fd, alignment_file = tempfile.mkstemp(suffix='.sto')
  os.close(fd)
  try:
    subprocess.check_call(
      'cmalign -g --notrunc {} {} > {}'.format(cm_model, ref_fasta, alignment_file),
      shell=True
    )
    ss, positions = parse_alignment_cons(alignment_file)
    col_to_sprinzl = build_col_to_sprinzl(ss, positions)
    sprinzl_axis, ref_to_sprinzl = build_sprinzl_mapping(alignment_file, col_to_sprinzl)
  finally:
    os.remove(alignment_file)

  return sprinzl_axis, ref_to_sprinzl


# ---------------------------------------------------------------------------
# Canonical Sprinzl order (used to reconstruct sprinzl_axis from a mapping)
# ---------------------------------------------------------------------------

# Labels in the order they appear in the consensus alignment columns.
# Insertions (e.g. '36i1') are not listed here; they sort immediately after
# their base position via build_axis_from_mapping().
_CANONICAL_SPRINZL = [
  '1', '2', '3', '4', '5', '6', '7', '8', '9',
  '10', '11', '12', '13', '14', '15', '16', '17', '17a', '18', '19',
  '20', '20a', '20b', '21', '22', '23', '24', '25', '26',
  '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '38',
  '39', '40', '41', '42', '43', '44', '45',
  'V11', 'V12', 'V13', 'V14', 'V15', 'V16', 'V17',
  'V1', 'V2', 'V3', 'V4', 'V5',
  'V27', 'V26', 'V25', 'V24', 'V23', 'V22', 'V21',
  '46', '47', '48', '49', '50', '51', '52', '53', '54', '55', '56',
  '57', '58', '59', '60', '61', '62', '63', '64', '65',
  '66', '67', '68', '69', '70', '71', '72', '73', '74', '75', '76',
]

_CANONICAL_RANK = {lbl: i for i, lbl in enumerate(_CANONICAL_SPRINZL)}


def _sprinzl_sort_key(label):
  """Sort key that places insertions (e.g. '36i1') after their base position."""
  if 'i' in label:
    base, ins = label.rsplit('i', 1)
    return (_CANONICAL_RANK.get(base, len(_CANONICAL_SPRINZL)), 1, int(ins))
  return (_CANONICAL_RANK.get(label, len(_CANONICAL_SPRINZL)), 0, 0)


def build_axis_from_mapping(ref_to_sprinzl):
  """
  Derive an ordered sprinzl_axis from the union of labels in ref_to_sprinzl.

  Used when loading a hand-edited mapping TSV without running cmalign —
  the axis order is reconstructed from the canonical Sprinzl numbering rather
  than stored in the file.
  """
  all_labels = {lbl for labels in ref_to_sprinzl.values() for lbl in labels}
  return sorted(all_labels, key=_sprinzl_sort_key)


# ---------------------------------------------------------------------------
# Human-readable mapping I/O
# ---------------------------------------------------------------------------

def save_sprinzl_mapping(sprinzl_axis, ref_to_sprinzl, path):
  """
  Write ref_to_sprinzl to a tab-separated file for inspection or hand-editing.

  Format: header row (ref_name / ref_position / sprinzl_label), then one row
  per physical reference position (1-based ref_position).  Opens directly in
  Excel/LibreOffice with no import options needed.  Reload with
  load_sprinzl_mapping() via --sprinzl-map; axis order is reconstructed from
  the canonical Sprinzl numbering automatically.
  """
  with open(path, 'w') as fh:
    fh.write('ref_name\tref_position\tsprinzl_label\n')
    for ref_name, labels in ref_to_sprinzl.items():
      for pos_idx, label in enumerate(labels, start=1):
        fh.write('{}\t{}\t{}\n'.format(ref_name, pos_idx, label))


def load_sprinzl_mapping(path):
  """
  Load a Sprinzl mapping TSV written by save_sprinzl_mapping (or hand-edited).

  Parameters
  ----------
  path : str
      Path to the TSV file.

  Returns
  -------
  ref_to_sprinzl : dict[str, list[str]]
      {ref_name: [sprinzl_label, ...]} sorted by ref_position (1-based).
  """
  rows = {}  # {ref_name: {ref_position(int): sprinzl_label}}

  with open(path) as fh:
    for line in fh:
      line = line.rstrip('\n')
      if line.startswith('#') or not line:
        continue
      parts = line.split('\t')
      if len(parts) != 3:
        continue
      ref_name, ref_pos_str, label = parts
      if ref_name == 'ref_name':  # skip header row
        continue
      try:
        ref_pos = int(ref_pos_str)
      except ValueError:
        continue
      if ref_name not in rows:
        rows[ref_name] = {}
      rows[ref_name][ref_pos] = label

  return {
    ref_name: [label for _, label in sorted(pos_map.items())]
    for ref_name, pos_map in rows.items()
  }
