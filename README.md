# tRNA-heatmap

Visualizes tRNA alignment pileup data as a mismatch-rate heatmap with a shared Sprinzl coordinate x-axis.

---

## Requirements

- Python >= 3.9
- Python packages: `numpy`, `pysam`, `pandas`, `matplotlib`
- External tool: [Infernal](http://eddylab.org/infernal/) (`cmalign`) — required for the `run` and `inspect` subcommands

---

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/genometechlab/tRNA-heatmap/
   cd tRNA-heatmap
   ```

2. Install the Python package and its dependencies:
   ```bash
   pip install -e .
   ```

3. Install Infernal (required for `run` and `inspect`):
   ```bash
   conda install -c bioconda infernal
   # or, on macOS with Homebrew:
   brew install infernal
   ```

After installation, the `tRNAheatmap` command will be available in your terminal.

---

## Quick Start

The `run` subcommand is the full pipeline. All input BAMs are organized into named `--condition` groups: pass one group to get a single heatmap, two or more to get pairwise delta heatmaps in one step.

**Single heatmap (one condition)**

```bash
# Using a bundled organism model (recommended)
tRNAheatmap run \
  --ref ref.fa \
  --condition sample aligned.bam \
  --organism eukaryotic \
  --output heatmap.png \
  --save-df pileup \
  --threads 4

# Using a custom covariance model
tRNAheatmap run \
  --ref ref.fa \
  --condition sample aligned.bam \
  --cm /path/to/my-model.cm \
  --output heatmap.png \
  --save-df pileup \
  --threads 4
```

**Pairwise delta heatmaps (≥2 conditions)**

```bash
# Two single-BAM conditions
tRNAheatmap run \
  --ref ref.fa \
  --condition condA condA.bam \
  --condition condB condB.bam \
  --organism eukaryotic \
  --output results/delta \
  --threads 4

# Replicate-aware: average per-BAM mismatch rates within each condition
tRNAheatmap run \
  --ref ref.fa \
  --organism eukaryotic --merge-mode equal \
  --condition WT wt1.bam wt2.bam wt3.bam \
  --condition KO ko1.bam ko2.bam ko3.bam \
  --output results/delta \
  --individual \
  --threads 4
```

With ≥2 `--condition` groups, `--output` is treated as a prefix: each pair is written as `{prefix}_{A}_vs_{B}.{ext}` (default extension `.pdf`). Adding `--individual` also emits one heatmap per condition (`{prefix}_{condition}.{ext}`).

| Flag | Description |
|------|-------------|
| `--ref` / `-r` | Reference tRNA FASTA used for alignment |
| `--condition` / `-C` | Condition group: first token is the name, remaining tokens are BAM paths. Repeat for multiple conditions. At least one is required. |
| `--merge-mode` | `total` (default): sum raw counts within a condition. `equal`: average per-BAM mismatch rates so each replicate weighs equally regardless of depth. |
| `--organism` / `-g` | Use a bundled model: `eukaryotic`, `prokaryotic`, or `archaeal` (mutually exclusive with `--cm`) |
| `--cm` / `-c` | Path to a custom Infernal covariance model `.cm` file (mutually exclusive with `--organism`) |
| `--output` / `-o` | One condition: output heatmap file (extension determines format). Multiple conditions: prefix; pairs become `{prefix}_{A}_vs_{B}.{ext}`. |
| `--save-df` / `-s` | Base name to save pileup data as a `.tsv` file. Only valid with exactly one `--condition`. |
| `--individual` | With ≥2 `--condition` groups, also emit a per-condition heatmap (`{prefix}_{condition}.{ext}`) alongside the pairwise deltas. |
| `--threads` / `-t` | Number of CPU threads (default: 1) |
| `--save-mapping BASE` | Save the computed Sprinzl coordinate mapping to `BASE.tsv` for inspection or hand-editing |
| `--sprinzl-map FILE` | Load a pre-computed or hand-edited mapping TSV; patches cmalign output if `--organism`/`--cm` is also given, or replaces it entirely if used alone |

**`inspect` — visualize Sprinzl coordinate coverage (no BAM needed)**

```bash
tRNAheatmap inspect \
  --ref ref.fa \
  --organism eukaryotic \
  --output sprinzl_coverage.png
```

Runs `cmalign` on the reference FASTA and renders a presence/absence heatmap: solid color = tRNA has a base at that Sprinzl position; black dot = no base (cmalign gap). No reads or pileup required. Use this to audit how Infernal handles a new reference set before committing to full runs.

| Flag | Description |
|------|-------------|
| `--ref` / `-r` | Reference tRNA FASTA file |
| `--organism` / `-g` | Bundled model: `eukaryotic`, `prokaryotic`, or `archaeal` (mutually exclusive with `--cm`) |
| `--cm` / `-c` | Path to a custom `.cm` file (mutually exclusive with `--organism`) |
| `--output` / `-o` | Output file; format set by extension. Default: `sprinzl_coverage.pdf` |
| `--save-mapping BASE` | Save the computed Sprinzl coordinate mapping to `BASE.tsv` for inspection or hand-editing |
| `--sprinzl-map FILE` | Load a pre-computed or hand-edited mapping TSV; patches cmalign output if `--organism`/`--cm` is also given, or replaces it entirely if used alone |

**Adapter / trimming flags (`run`, `inspect`)**

Use these to remove constant adapter sequences from the reference FASTA before Sprinzl coordinate assignment. The pileup arrays are trimmed to match.

| Flag | Description |
|------|-------------|
| `--detect-adapters` | Auto-detect and trim common 5′/3′ prefix/suffix shared by all tRNAs. Requires ≥2 sequences. Mutually exclusive with `--trim-5`. |
| `--trim-5 N` | Trim exactly N bases from the 5′ end of each reference. Mutually exclusive with `--detect-adapters`. |
| `--trim-3 N` | Trim exactly N bases from the 3′ end of each reference. Can be combined with `--trim-5`. |

**Reference filtering flags (`run`, `inspect`)**

| Flag | Description |
|------|-------------|
| `--include-refs FILE` | Plain-text file (one reference name per line). Only the listed references appear in the output. Mutually exclusive with `--exclude-refs`. |
| `--exclude-refs FILE` | Plain-text file (one reference name per line). Listed references are dropped. Mutually exclusive with `--include-refs`. |

**Shared plot options (all subcommands)**

| Flag | Default | Description |
|------|---------|-------------|
| `--palette` | `light-high` | Color scheme: `light-high` (yellow=low), `dark-high`, or `viridis` |
| `--title` | *(subcommand-specific)* | Figure title |
| `--ylabel` | `tRNA Reference Names` | Y-axis label |
| `--include-insertions` | off | Show insertion columns (e.g. `36i1`) in the plot |
| `--dpi` | 300 | Output resolution in DPI |
| `--cell-size` | 0.25 | Inches per heatmap cell for auto figure sizing |
| `--style FILE` | — | Matplotlib style sheet (`.mplstyle`) layered on top of the bundled default |
| `--outdir` / `-O` | cwd | Directory for all output files (created if missing) |

---

## Test Data Examples

The `test_data/` directory contains a yeast reference FASTA and 5%-subsampled BAM files for testing.

### 1. Single heatmap

Run the full pipeline on one BAM and save the pileup as a TSV for inspection.

```bash
tRNAheatmap run \
  --ref test_data/sacCer3-mature-tRNAs_zap_ref.fa \
  --condition WT test_data/yeast_wt_example_subsample_rep1.bam \
  --organism eukaryotic \
  --output results/wt_heatmap.png \
  --save-df results/wt \
  --threads 4
```

### 2. Inspect Sprinzl coordinate coverage (no BAM needed)

Audit how Infernal maps the reference sequences before running a full pileup.

```bash
tRNAheatmap inspect \
  --ref test_data/sacCer3-mature-tRNAs_zap_ref.fa \
  --organism eukaryotic \
  --output results/sprinzl_coverage.png
```

To detect and trim common adapter sequences automatically:

```bash
tRNAheatmap inspect \
  --ref test_data/sacCer3-mature-tRNAs_zap_ref.fa \
  --organism eukaryotic \
  --detect-adapters \
  --output results/sprinzl_coverage_trimmed.png
```

### 3. Delta heatmap between conditions (with replicates)

```bash
tRNAheatmap run \
  --ref test_data/sacCer3-mature-tRNAs_zap_ref.fa \
  --organism eukaryotic --merge-mode equal \
  --condition WT \
    test_data/yeast_wt_example_subsample_rep1.bam \
    test_data/yeast_wt_example_subsample_rep2.bam \
    test_data/yeast_wt_example_subsample_rep3.bam \
  --condition MT \
    test_data/yeast_mt_enrichment_example_subsample_rep1.bam \
    test_data/yeast_mt_enrichment_example_subsample_rep2.bam \
    test_data/yeast_mt_enrichment_example_subsample_rep3.bam \
  --output results/delta \
  --individual \
  --threads 4
```

This produces `results/delta_WT_vs_MT.pdf` (the pairwise delta) and, with `--individual`, `results/delta_WT.pdf` and `results/delta_MT.pdf` (the per-condition merged heatmaps).

---

## Overriding Sprinzl Coordinate Assignments

Infernal's `cmalign` assigns Sprinzl positions automatically, but unusual tRNA variants — extra loops, truncated stems, non-canonical anticodon regions — can receive incorrect assignments. When that happens, a domain expert can export the mapping, fix the affected rows, and feed the corrected file back in without rerunning the full pipeline.

### Step 1 — Export the mapping

Add `--save-mapping` to any `run` or `inspect` call. The flag writes `BASE.tsv` alongside the heatmap output.

```bash
tRNAheatmap run \
  --ref ref.fa \
  --condition sample aligned.bam \
  --organism eukaryotic \
  --output heatmap.png \
  --save-mapping my_coords
```

This produces `my_coords.tsv`. If you only want to inspect the mapping without running a pileup, use `inspect` instead — it is faster because it skips BAM reading entirely:

```bash
tRNAheatmap inspect \
  --ref ref.fa \
  --organism eukaryotic \
  --output sprinzl_coverage.png \
  --save-mapping my_coords
```

### Step 2 — Inspect and edit the TSV

The file looks like this:

```
ref_name	ref_position	sprinzl_label	modification_symbol
tRNA-Ala-AGC-1	1	1	
tRNA-Ala-AGC-1	2	2	
...
tRNA-Ala-AGC-1	34	34	
tRNA-Ala-AGC-1	35	35	
...
```

The columns are tab-separated: sequence name, reference position (1-based), Sprinzl label, and an optional modification symbol. To correct a mis-assigned position, change only the `sprinzl_label` column. For example, if position 34 of `tRNA-Ala-AGC-1` was assigned label `34` but should be `33`:

```
tRNA-Ala-AGC-1	34	33	
```

Leave all other rows unchanged. The `modification_symbol` column is described below.

### Step 3 — Re-run with the edited mapping

**Patch mode** — use the edited TSV alongside `--organism` or `--cm`. Only the sequences listed in the TSV are overwritten; all other references are computed normally by cmalign.

```bash
tRNAheatmap run \
  --ref ref.fa \
  --condition sample aligned.bam \
  --organism eukaryotic \
  --sprinzl-map my_coords.tsv \
  --output heatmap_corrected.png
```

**Full replacement mode** — omit `--organism` and `--cm`. cmalign is skipped entirely. The TSV must cover every reference sequence in the FASTA.

```bash
tRNAheatmap run \
  --ref ref.fa \
  --condition sample aligned.bam \
  --sprinzl-map my_coords.tsv \
  --output heatmap_corrected.png
```

---

## Annotating Modified Nucleotides

The mapping TSV supports an optional fourth column, `modification_symbol`, for marking known RNA modifications. When present, the symbol (any Unicode character — e.g. `ψ` for pseudouridine, `I` for inosine, `m` for m⁶A) is drawn centered on the corresponding heatmap cell. The text color is chosen automatically for contrast against the cell's background color.

To annotate modifications, add the symbol in the `modification_symbol` column for the relevant rows and re-run with `--sprinzl-map`:

```
ref_name	ref_position	sprinzl_label	modification_symbol
tRNA-Ala-AGC-1	34	34	ψ
tRNA-Ala-AGC-1	37	37	I
```

Rows with no modification can leave the fourth column empty or omit it entirely. The modification symbol is rendered on top of the mismatch-rate color; cells with no coverage (NaN) or a cmalign gap (black dot) are skipped.

---

## Input File Requirements

You need three inputs for the `run` subcommand: a **reference FASTA** containing your tRNA sequences, one or more **sorted BAM files** of reads aligned to that FASTA (produced by any short- or long-read aligner) grouped into named `--condition` arguments, and a **covariance model**. Use `--organism eukaryotic/prokaryotic/archaeal` to select a bundled model automatically, or `--cm /path/to/model.cm` to supply your own. The bundled models are installed with the package — no path resolution needed.

> **Important:** The sequence names in your BAM files and the names used by `cmalign` must match — both must come from the same reference FASTA. If you see missing data in the heatmap, verify name consistency with:
> ```bash
> samtools view -H aligned.bam | grep '^@SQ'
> ```

---

## Output

| File | Description |
|------|-------------|
| `heatmap.png` (or `.pdf`/`.svg`) | Heatmap with rows = tRNA references, columns = Sprinzl positions. Color = mismatch rate (0–1). Blank cell = position exists but had zero reads. Black dot (●) = cmalign gap (no base at that Sprinzl position). |
| `{prefix}_{A}_vs_{B}.{ext}` | Pairwise delta heatmap (multi-condition runs). Diverging Spectral_r colormap, range −1…+1. Black dot = either condition lacks the position. |
| `{prefix}_{condition}.{ext}` | Per-condition merged heatmap (multi-condition runs with `--individual`). |
| `pileup.tsv` (from `--save-df` in `total` mode) | Long-format table with per-position `match`, `mismatch`, `insertion`, `deletion`, `accuracy`, and `mismatch_rate` columns. |
| `rates.tsv` (from `--save-df` in `equal` mode) | Long-format table with `mismatch_rate` and `std_dev` columns (std is NaN for single-BAM conditions). |
| `coords.tsv` (from `--save-mapping`) | Sprinzl coordinate mapping exported by Infernal. Tab-separated columns: `ref_name`, `ref_position` (1-based), `sprinzl_label`, and optional `modification_symbol`. Opens directly in Excel or LibreOffice and can be hand-edited before reuse with `--sprinzl-map`. |

---

## Acknowledgements

The `.cm` covariance model files and Sprinzl position mapping code are adapted from [tRNAviz](https://github.com/UCSC-LoweLab/tRNAviz) (LGPL-3.0) by the UCSC Lowe Lab (Lin BY, Chan PP, Lowe TM. 2019. tRNAviz. *Nucleic Acids Res.*).

This tool was developed with assistance from [Claude Code](https://claude.ai/code) (Anthropic).
