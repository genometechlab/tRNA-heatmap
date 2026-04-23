# tRNA-heatmap

Visualizes tRNA alignment pileup data as a mismatch-rate heatmap with a shared Sprinzl coordinate x-axis.

---

## Requirements

- Python >= 3.9
- Python packages: `numpy`, `pysam`, `pandas`, `matplotlib`
- External tool: [Infernal](http://eddylab.org/infernal/) (`cmalign`) — required only for the `run` subcommand

---

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/your-org/tRNA-heatmap.git
   cd tRNA-heatmap
   ```

2. Install the Python package and its dependencies:
   ```bash
   pip install -e .
   ```

3. Install Infernal (skip if you only plan to use `plot` or `delta` with pre-computed data):
   ```bash
   conda install -c bioconda infernal
   # or, on macOS with Homebrew:
   brew install infernal
   ```

After installation, the `tRNAheatmap` command will be available in your terminal.

---

## Quick Start

**`run` — full pipeline (BAM to heatmap)**

```bash
# Using a bundled organism model (recommended)
tRNAheatmap run \
  --ref ref.fa \
  --bam aligned.bam \
  --organism eukaryotic \
  --output heatmap.png \
  --save-df pileup \
  --threads 4

# Using a custom covariance model
tRNAheatmap run \
  --ref ref.fa \
  --bam aligned.bam \
  --cm /path/to/my-model.cm \
  --output heatmap.png \
  --save-df pileup \
  --threads 4
```

| Flag | Description |
|------|-------------|
| `--ref` / `-r` | Reference tRNA FASTA used for alignment |
| `--bam` / `-b` | Sorted BAM file of reads aligned to the reference |
| `--organism` / `-g` | Use a bundled model: `eukaryotic`, `prokaryotic`, or `archaeal` (mutually exclusive with `--cm`) |
| `--cm` / `-c` | Path to a custom Infernal covariance model `.cm` file (mutually exclusive with `--organism`) |
| `--output` / `-o` | Output file; format set by extension (`.png`, `.pdf`, `.svg`) |
| `--save-df` / `-s` | Base name to save pileup data as `.npz` (for replotting) and `.tsv` (for inspection) |
| `--threads` / `-t` | Number of CPU threads (default: 1) |
| `--save-mapping BASE` | Save the computed Sprinzl coordinate mapping to `BASE.tsv` for inspection or hand-editing |
| `--sprinzl-map FILE` | Load a pre-computed or hand-edited mapping TSV; patches cmalign output if `--organism`/`--cm` is also given, or replaces it entirely if used alone |

**`plot` — replot from saved data (no BAM or Infernal needed)**

```bash
tRNAheatmap plot \
  --load pileup.npz \
  --output heatmap.png
```

**`run-delta` — full pipeline straight to a delta heatmap**

```bash
tRNAheatmap run-delta \
  --ref ref.fa \
  --bam condA.bam condB.bam \
  --organism eukaryotic \
  --output results/delta \
  --threads 4
```

Runs pileup on each BAM, computes Sprinzl coordinates once, and produces one delta heatmap per pair — no intermediate NPZ files required. Use `--bam` with three or more files to compare all pairs at once.

| Flag | Description |
|------|-------------|
| `--ref` / `-r` | Shared reference tRNA FASTA |
| `--bam` / `-b` | Two or more sorted BAM files to compare |
| `--organism` / `-g` | Bundled model: `eukaryotic`, `prokaryotic`, or `archaeal` (mutually exclusive with `--cm`) |
| `--cm` / `-c` | Path to a custom `.cm` file (mutually exclusive with `--organism`) |
| `--output` / `-o` | Output prefix; produces `{prefix}_{A}_vs_{B}.{ext}` per pair. Default: `delta` |
| `--threads` / `-t` | Number of CPU threads (default: 1) |
| `--save-mapping BASE` | Save the computed Sprinzl coordinate mapping to `BASE.tsv` for inspection or hand-editing |
| `--sprinzl-map FILE` | Load a pre-computed or hand-edited mapping TSV; patches cmalign output if `--organism`/`--cm` is also given, or replaces it entirely if used alone |

**`delta` — compare two or more pre-computed samples**

```bash
tRNAheatmap delta \
  --load sampleA.npz sampleB.npz \
  --output results/delta
```

Produces one difference heatmap per pair (mismatch rate A minus mismatch rate B), saved as `delta_sampleA_vs_sampleB.png`. Use this when you already have NPZ files from `run --save-df`.

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
| `--output` / `-o` | Output file; format set by extension. Default: `sprinzl_coverage.png` |
| `--save-mapping BASE` | Save the computed Sprinzl coordinate mapping to `BASE.tsv` for inspection or hand-editing |
| `--sprinzl-map FILE` | Load a pre-computed or hand-edited mapping TSV; patches cmalign output if `--organism`/`--cm` is also given, or replaces it entirely if used alone |

**Adapter / trimming flags (`run`, `run-delta`, `inspect`)**

Use these to remove constant adapter sequences from the reference FASTA before Sprinzl coordinate assignment. The pileup arrays are trimmed to match.

| Flag | Description |
|------|-------------|
| `--detect-adapters` | Auto-detect and trim common 5′/3′ prefix/suffix shared by all tRNAs. Requires ≥2 sequences. Mutually exclusive with `--trim-5`. |
| `--trim-5 N` | Trim exactly N bases from the 5′ end of each reference. Mutually exclusive with `--detect-adapters`. |
| `--trim-3 N` | Trim exactly N bases from the 3′ end of each reference. Can be combined with `--trim-5`. |

**Shared plot options (all subcommands)**

| Flag | Default | Description |
|------|---------|-------------|
| `--palette` | `light-high` | Color scheme: `light-high` (yellow=low), `dark-high`, or `viridis` |
| `--title` | *(subcommand-specific)* | Figure title |
| `--ylabel` | `tRNA Reference Names` | Y-axis label |
| `--include-insertions` | off | Show insertion columns (e.g. `36i1`) in the plot |
| `--dpi` | 300 | Output resolution in DPI |

---

## Test Data Examples

The `test_data/` directory contains a yeast reference FASTA and two 5%-subsampled BAM files (condition A: deacylated; condition B: deacylless CCA) for testing.

### 1. Single heatmap

Run the full pipeline on condition A and save the pileup for later replotting.

```bash
tRNAheatmap run \
  --ref test_data/sacCer3-mature-tRNAs_zap_ref.fa \
  --bam test_data/09_11_25_RNA004_Bio_tRNA_deAcy_Yeast_mt.dorado.1.0.0.emit_moves.sorted.5percent_subsample.bam \
  --organism eukaryotic \
  --output results/condA_heatmap.png \
  --save-df results/condA \
  --threads 4
```

Replot from the saved NPZ without rerunning the pileup or `cmalign`:

```bash
tRNAheatmap plot \
  --load results/condA.npz \
  --output results/condA_heatmap_replot.png
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

### 3. Delta heatmap between two conditions

**One-step approach using `run-delta`:**

```bash
tRNAheatmap run-delta \
  --ref test_data/sacCer3-mature-tRNAs_zap_ref.fa \
  --bam test_data/09_11_25_RNA004_Bio_tRNA_deAcy_Yeast_mt.dorado.1.0.0.emit_moves.sorted.5percent_subsample.bam \
       test_data/09_15_25_sctRNAdeaclessCCA.dorado_1.0.0.emit_moves.sorted.5percent_subsample.bam \
  --organism eukaryotic \
  --output results/delta \
  --threads 4
```

**Three-step approach using `run` + `delta`** (use this when you want to reuse the NPZ files later):

```bash
tRNAheatmap run \
  --ref test_data/sacCer3-mature-tRNAs_zap_ref.fa \
  --bam test_data/09_11_25_RNA004_Bio_tRNA_deAcy_Yeast_mt.dorado.1.0.0.emit_moves.sorted.5percent_subsample.bam \
  --organism eukaryotic \
  --output results/condA_heatmap.png \
  --save-df results/condA \
  --threads 4

tRNAheatmap run \
  --ref test_data/sacCer3-mature-tRNAs_zap_ref.fa \
  --bam test_data/09_15_25_sctRNAdeaclessCCA.dorado_1.0.0.emit_moves.sorted.5percent_subsample.bam \
  --organism eukaryotic \
  --output results/condB_heatmap.png \
  --save-df results/condB \
  --threads 4

tRNAheatmap delta \
  --load results/condA.npz results/condB.npz \
  --output results/delta
```

---

## Overriding Sprinzl Coordinate Assignments

Infernal's `cmalign` assigns Sprinzl positions automatically, but unusual tRNA variants — extra loops, truncated stems, non-canonical anticodon regions — can receive incorrect assignments. When that happens, a domain expert can export the mapping, fix the affected rows, and feed the corrected file back in without rerunning the full pipeline.

### Step 1 — Export the mapping

Add `--save-mapping` to any `run`, `run-delta`, or `inspect` call. The flag writes `BASE.tsv` alongside the heatmap output.

```bash
tRNAheatmap run \
  --ref ref.fa \
  --bam aligned.bam \
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
#axis:1,2,3,...,34,35,36,...,73,74,75
tRNA-Ala-AGC-1	1	1
tRNA-Ala-AGC-1	2	2
...
tRNA-Ala-AGC-1	34	34
tRNA-Ala-AGC-1	35	35
...
```

The first line (`#axis:`) records the ordered global set of Sprinzl labels — all subcommands need this to align the shared x-axis. The remaining lines are tab-separated: sequence name, reference position (1-based), Sprinzl label.

To correct a mis-assigned position, change only the `sprinzl_label` column. For example, if position 34 of `tRNA-Ala-AGC-1` was assigned label `34` but should be `33`:

```
tRNA-Ala-AGC-1	34	33
```

Leave all other rows and the `#axis:` header unchanged.

### Step 3 — Re-run with the edited mapping

**Patch mode** — use the edited TSV alongside `--organism` or `--cm`. Only the sequences listed in the TSV are overwritten; all other references are computed normally by cmalign.

```bash
tRNAheatmap run \
  --ref ref.fa \
  --bam aligned.bam \
  --organism eukaryotic \
  --sprinzl-map my_coords.tsv \
  --output heatmap_corrected.png
```

**Full replacement mode** — omit `--organism` and `--cm`. cmalign is skipped entirely. The TSV must include the `#axis:` header and must cover every reference sequence in the FASTA.

```bash
tRNAheatmap run \
  --ref ref.fa \
  --bam aligned.bam \
  --sprinzl-map my_coords.tsv \
  --output heatmap_corrected.png
```

---

## Input File Requirements

You need three inputs for the `run` subcommand: a **reference FASTA** containing your tRNA sequences, a **sorted BAM file** of reads aligned to that FASTA (produced by any short- or long-read aligner), and a **covariance model**. Use `--organism eukaryotic/prokaryotic/archaeal` to select a bundled model automatically, or `--cm /path/to/model.cm` to supply your own. The bundled models are installed with the package — no path resolution needed.

> **Important:** The sequence names in your BAM file and the names used by `cmalign` must match — both must come from the same reference FASTA. If you see missing data in the heatmap, verify name consistency with:
> ```bash
> samtools view -H aligned.bam | grep '^@SQ'
> ```

---

## Output

| File | Description |
|------|-------------|
| `heatmap.png` (or `.pdf`/`.svg`) | Heatmap with rows = tRNA references, columns = Sprinzl positions. Color = mismatch rate (0–1). Blank cell = position exists but had zero reads. Black dot (●) = cmalign gap (no base at that Sprinzl position). |
| `pileup.npz` (from `--save-df`) | Saved pileup array for fast replotting with `plot` or `delta`. |
| `pileup.tsv` (from `--save-df`) | Long-format table with per-position `match`, `mismatch`, `insertion`, `deletion`, `accuracy`, and `mismatch_rate` columns. |
| `coords.tsv` (from `--save-mapping`) | Sprinzl coordinate mapping exported by Infernal. Contains a `#axis:` comment header listing the global Sprinzl label order, followed by tab-separated rows of `ref_name`, `ref_position` (1-based), and `sprinzl_label`. Opens directly in Excel or LibreOffice and can be hand-edited before reuse with `--sprinzl-map`. |

---

## Acknowledgements

The `.cm` covariance model files and Sprinzl position mapping code are adapted from [tRNAviz](https://github.com/UCSC-LoweLab/tRNAviz) (LGPL-3.0) by the UCSC Lowe Lab (Lin BY, Chan PP, Lowe TM. 2019. tRNAviz. *Nucleic Acids Res.*).

This tool was developed with assistance from [Claude Code](https://claude.ai/code) (Anthropic).
