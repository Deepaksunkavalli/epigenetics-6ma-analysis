# Epigenetics 6mA Analysis

This repository contains scripts and derived data for a small‑scale exploration of N6‑methyladenine (6mA) motifs in the *Arabidopsis thaliana* genome.

## Contents

- `6ma_motif_analysis.py` – self‑contained Python pipeline that scans a genome FASTA for a user‑defined motif (default: GAGG), extracts flanking windows, summarizes motif counts per chromosome, computes nucleotide frequency matrices and produces a sequence logo.
- `report.md` – research report summarizing findings from scanning the Col‑CEN assembly and integrating histone modification datasets.
- `motif_counts_by_chromosome.tsv` – tab‑delimited summary of motif counts per chromosome for the Col‑CEN assembly.
- `motif_position_logo.png` – sequence logo plot for the GAGG motif.
- `GAGG_motif_positions.bed` – BED file of all GAGG motif occurrences in the Col‑CEN assembly.
- `GAGG_windows.fasta.gz` – compressed FASTA sequences (±30 bp) around each GAGG motif occurrence.
- `motif_position_frequency.npy` – NumPy array containing the nucleotide frequency matrix used to build the sequence logo.
- `dataset_sources.md` – documentation describing where to obtain the raw genome assemblies and histone datasets.

## Usage

To run the motif analysis on another genome assembly or motif:

```
python 6ma_motif_analysis.py --fasta path/to/genome.fasta --motif MOTIF --flank 30 --outdir results
```

Replace `MOTIF` with the motif of interest and adjust `--flank` as needed. The script will scan the genome for the motif, extract windows of length `2×flank + len(MOTIF)`, and write output files to the specified directory.

## Data sources

The raw genome assemblies (TAIR10, Col‑CEN) and histone modification peak sets are not included in this repository due to their size. Please refer to `dataset_sources.md` for links and instructions on obtaining these datasets.
