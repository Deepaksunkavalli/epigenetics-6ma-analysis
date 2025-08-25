#!/usr/bin/env python3
"""
6mA motif analysis pipeline for Arabidopsis thaliana genomes.

This script implements a self‑contained analysis pipeline for scanning a genome
FASTA file for the consensus N6‑methyladenine (6mA) motif "GAGG" and its
reverse complement "CCTC".  It extracts genomic windows around each
occurrence, summarises motif distributions and generates output files ready
for downstream epigenetic analyses.

The pipeline is designed to work without external dependencies beyond
standard Python packages (pandas, numpy, matplotlib).  It does not rely
on Biopython or other bioinformatics libraries, making it suitable for
restricted environments.

Key output files:
  - GAGG_motif_positions.bed: BED format intervals for all motif sites.
  - GAGG_windows.fasta: FASTA sequences surrounding motif sites.
  - motif_counts_by_chromosome.tsv: Tab‑delimited summary of motif counts per chromosome.
  - motif_position_frequency.npy: NumPy array of nucleotide frequencies across window positions.
  - motif_position_logo.png: Sequence logo visualising base frequencies across the motif window.

Usage:
  python 6ma_motif_analysis.py --fasta Col-CEN_v1.2.fasta --window 15 --outdir results

This will scan the provided genome, extract 30 bp windows (15 bp upstream and
downstream around the 4 bp motif) and write all outputs to the specified
directory.

Limitations:
  * The script scans a single genome FASTA file.  For comparative analyses
    between assemblies (e.g., TAIR10 vs. Col‑PEK), run the script separately
    for each assembly and compare the resulting BED or TSV files using
    BEDTools or pandas.
  * Histone modification integration is not performed here but can be
    accomplished by intersecting the BED file with histone peak BED files.

"""

import argparse
import os
from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def parse_fasta(fasta_path):
    """Generator that yields (header, sequence) tuples from a FASTA file.

    Parameters
    ----------
    fasta_path : str
        Path to the FASTA file to parse.

    Yields
    ------
    tuple
        (header, sequence) where header is the line after '>' and sequence is a
        concatenated upper‑case string containing only ACGT (other bases are
        retained as given).
    """
    header = None
    seq_chunks = []
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                # Emit previous record
                if header is not None:
                    yield header, ''.join(seq_chunks).upper()
                    seq_chunks = []
                header = line[1:]
            else:
                seq_chunks.append(line)
        # Emit the last record
        if header is not None:
            yield header, ''.join(seq_chunks).upper()


def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    complement = str.maketrans('ACGTN', 'TGCAN')
    return seq.translate(complement)[::-1]


def scan_motif(chrom, sequence, motif='GAGG'):
    """Scan a chromosome sequence for occurrences of a motif and its reverse complement.

    Parameters
    ----------
    chrom : str
        Chromosome name.
    sequence : str
        Upper‑case DNA sequence.
    motif : str, optional
        6mA consensus motif to scan (forward strand). The reverse complement
        will automatically be derived. Default is 'GAGG'.

    Returns
    -------
    list of dict
        Each dict has keys: 'chrom', 'start', 'end', 'strand'. Coordinates
        are zero‑based [start, end) intervals following BED convention.
    """
    motif = motif.upper()
    rc_motif = reverse_complement(motif)
    motif_len = len(motif)
    results = []
    # Scan forward strand
    idx = sequence.find(motif)
    while idx != -1:
        results.append({'chrom': chrom, 'start': idx, 'end': idx + motif_len, 'strand': '+'})
        idx = sequence.find(motif, idx + 1)
    # Scan reverse complement
    if rc_motif != motif:  # avoid double counting palindromic motifs
        idx = sequence.find(rc_motif)
        while idx != -1:
            results.append({'chrom': chrom, 'start': idx, 'end': idx + motif_len, 'strand': '-'})
            idx = sequence.find(rc_motif, idx + 1)
    return results


def extract_windows(sequence, positions, flank):
    """Extract windows around motif positions.

    Parameters
    ----------
    sequence : str
        Full chromosome sequence.
    positions : list of dict
        List of motif occurrences with keys 'start', 'end', 'strand'.
    flank : int
        Number of bases upstream and downstream to include.

    Returns
    -------
    list of dict
        Each dict contains 'header' and 'seq' for a FASTA entry.
    """
    windows = []
    seq_len = len(sequence)
    for pos in positions:
        start = pos['start']
        end = pos['end']
        # define window boundaries (0‑based, inclusive for start, exclusive for end)
        win_start = max(0, start - flank)
        win_end = min(seq_len, end + flank)
        window_seq = sequence[win_start:win_end]
        # If motif on negative strand, take reverse complement of the window
        if pos['strand'] == '-':
            window_seq = reverse_complement(window_seq)
        header = f"{pos['chrom']}:{start}-{end}:{pos['strand']}"
        windows.append({'header': header, 'seq': window_seq})
    return windows


def write_bed(motif_positions, path):
    """Write motif positions to a BED file."""
    with open(path, 'w') as bed:
        for pos in motif_positions:
            bed.write(f"{pos['chrom']}	{pos['start']}	{pos['end']}	GAGG	0	{pos['strand']}\n")


def write_fasta(windows, path):
    """Write extracted windows to a FASTA file."""
    with open(path, 'w') as fasta:
        for entry in windows:
            fasta.write(f">{entry['header']}\n{entry['seq']}\n")


def compute_nucleotide_frequencies(windows, flank, motif_len):
    """Compute nucleotide frequency matrix for alignment windows.

    Parameters
    ----------
    windows : list of dict
        Extracted windows with 'seq' key. All sequences assumed to be same length.
    flank : int
        Number of flanking bases around motif.
    motif_len : int
        Length of the motif.

    Returns
    -------
    numpy.ndarray
        Array of shape (window_len, 4) containing frequencies of [A, C, G, T] at
        each position.
    """
    # Determine window length
    if not windows:
        return np.zeros((0, 4))
    window_len = len(windows[0]['seq'])
    counts = np.zeros((window_len, 4), dtype=int)
    nuc_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for entry in windows:
        seq = entry['seq']
        if len(seq) != window_len:
            # pad shorter sequences with Ns (ignored in counts)
            seq = seq.ljust(window_len, 'N')
        for i, base in enumerate(seq.upper()):
            if base in nuc_to_idx:
                counts[i, nuc_to_idx[base]] += 1
    frequencies = counts / counts.sum(axis=1, keepdims=True)
    return frequencies


def plot_sequence_logo(freq_matrix, path, motif_len, flank):
    """Create a simple sequence logo plot from a nucleotide frequency matrix.

    The plot shows stacked bars at each position representing the frequencies of
    A, C, G and T.  The motif region is highlighted with a shaded box.

    Parameters
    ----------
    freq_matrix : numpy.ndarray
        Array of shape (window_len, 4) with frequencies of A, C, G, T.
    path : str
        Output path for the PNG plot.
    motif_len : int
        Length of the motif (to delineate the motif region).
    flank : int
        Flanking length on each side of the motif.
    """
    # Determine number of positions
    window_len = freq_matrix.shape[0]
    x = np.arange(window_len)
    base_order = ['A', 'C', 'G', 'T']
    colors = {'A': '#2E7D32', 'C': '#0288D1', 'G': '#FBC02D', 'T': '#E64A19'}
    plt.figure(figsize=(12, 4))
    bottom = np.zeros(window_len)
    for i, base in enumerate(base_order):
        plt.bar(x, freq_matrix[:, i], bottom=bottom, color=colors[base], label=base)
        bottom += freq_matrix[:, i]
    # Highlight motif region
    motif_start = flank
    motif_end = flank + motif_len
    plt.axvspan(motif_start - 0.5, motif_end - 0.5, color='gray', alpha=0.2)
    plt.xlabel('Position relative to motif')
    plt.ylabel('Base frequency')
    plt.title('Nucleotide frequencies around 6mA motif GAGG')
    plt.legend(title='Base', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(path)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Scan genome for 6mA motif and generate summary files.")
    parser.add_argument('--fasta', required=True, help="Input genome FASTA file.")
    parser.add_argument('--motif', default='GAGG', help="Consensus motif to search (default: GAGG).")
    parser.add_argument('--flank', type=int, default=15, help="Number of bases upstream and downstream to include in windows.")
    parser.add_argument('--outdir', default='results', help="Directory to write output files.")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    all_positions = []
    windows = []
    motif_len = len(args.motif)

    # Scan each chromosome
    for header, seq in parse_fasta(args.fasta):
        # Use first token of header as chromosome name (strip anything after whitespace)
        chrom = header.split()[0]
        positions = scan_motif(chrom, seq, args.motif)
        all_positions.extend(positions)
        windows += extract_windows(seq, positions, args.flank)

    # Write BED file
    bed_path = os.path.join(args.outdir, 'GAGG_motif_positions.bed')
    write_bed(all_positions, bed_path)
    # Write FASTA of windows
    fasta_path = os.path.join(args.outdir, 'GAGG_windows.fasta')
    write_fasta(windows, fasta_path)
    # Summarise counts per chromosome
    df = pd.DataFrame(all_positions)
    counts = df.groupby('chrom').size().reset_index(name='count')
    counts_path = os.path.join(args.outdir, 'motif_counts_by_chromosome.tsv')
    counts.to_csv(counts_path, sep='\t', index=False)
    # Compute frequency matrix
    freq_matrix = compute_nucleotide_frequencies(windows, args.flank, motif_len)
    freq_path = os.path.join(args.outdir, 'motif_position_frequency.npy')
    np.save(freq_path, freq_matrix)
    # Plot sequence logo
    logo_path = os.path.join(args.outdir, 'motif_position_logo.png')
    if freq_matrix.size > 0:
        plot_sequence_logo(freq_matrix, logo_path, motif_len, args.flank)
    print(f"Processed {len(df)} motif occurrences across {counts.shape[0]} chromosomes.")
    print(f"Results saved in directory: {args.outdir}")


if __name__ == '__main__':
    main()