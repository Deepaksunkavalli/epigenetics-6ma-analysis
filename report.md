# Exploring N6‑methyladenine (6mA) motifs in *Arabidopsis thaliana*

## Background

N6‑methyladenine (6mA) is a DNA modification found in bacteria and some
eukaryotes.  In plants, the presence and function of 6mA has only recently
been appreciated.  Genome‑wide profiling revealed that 6mA marks are not
randomly distributed; instead, they occur preferentially at a short consensus
sequence.  A review on 6mA noted that plant genomes share a conserved
**GAGG** motif around methylated adenine and that this motif is also
observed in other organisms【733812400672236†L568-L575】.  The non‑random
placement of 6mA suggests a potential regulatory role in gene expression.

The reference genome for *Arabidopsis thaliana* has historically been the
TAIR10 assembly.  New assemblies such as Col‑XJTU and Col‑CEN (also referred
to as Col‑PEK/Col‑CEN) provide more contiguous sequences and additional
sequence compared with TAIR10.  The Col‑XJTU assembly spans **133 Mb** and
adds **~14.6 Mb** of novel sequence relative to TAIR10 while improving
centromeric and telomeric regions【825213852822405†L278-L297】.  Similarly,
the Col‑CEN assembly reconstructs the five centromeres, adds **12.6 Mb** of
new sequence and corrects presence/absence variation (27 genes missing
relative to TAIR10)【506180248870255†L598-L606】.  These assemblies enable
more accurate mapping of epigenetic marks.

Large collections of histone modification data are available for
*Arabidopsis*.  The ReMap project curated **795 ChIP‑seq/DAP‑seq datasets**,
including **286 histone datasets**, mapped to TAIR10; these datasets
represent **∼4.5 million** high‑quality histone peaks and are freely
available【404858133361378†L43-L97】.  These epigenomic tracks, together with
high‑quality genomes, provide a framework to explore relationships between
6mA motifs and chromatin states.

## Data sources

The Col‑CEN assembly (version 1.2) FASTA file was used as the primary
reference genome.  This assembly resolves centromeres and telomeres
and includes additional sequence compared with TAIR10【506180248870255†L598-L606】.
Although the TAIR10 genome is widely available via Ensembl and NCBI, network
restrictions prevented direct download, so the analysis focused on the
accessible Col‑CEN assembly.  The pipeline described below can be repeated
for TAIR10 or other assemblies by providing the appropriate FASTA file.

Histone modification data were not downloaded in this environment due to
similar restrictions.  Instead, instructions are provided for integrating
motif positions with histone peaks using BEDTools once the user obtains
the BED files (e.g. from the ReMap project【404858133361378†L43-L97】).

## Methods

The analysis was implemented in the script `6ma_motif_analysis.py` (see
GitHub repository).  The key steps are summarised below:

1. **Scanning for the 6mA consensus motif.**  The script reads each
   chromosome from a FASTA file and searches for occurrences of the
   consensus **GAGG** motif and its reverse complement **CCTC** (to
   capture motifs on either strand).  For every match, the chromosome,
   zero‑based start and end positions and strand are recorded.  The results
   are written to a BED file (`GAGG_motif_positions.bed`).

2. **Extracting sequence windows.**  For each motif occurrence, a window
   containing 15 bp upstream and 15 bp downstream around the 4‑bp motif
   (total 34 bp) is extracted.  For reverse‑strand hits, the window is
   reverse‑complemented to align motifs in the same orientation.  These
   windows are written to a FASTA file (`GAGG_windows.fasta`) for input
   to motif discovery tools such as MEME‑CHIP.

3. **Summary statistics.**  Motif counts are aggregated per chromosome,
   generating a simple table of motif density.  A nucleotide frequency
   matrix across the 34 bp windows is computed and saved as a NumPy
   array.  Using the frequencies, a sequence logo is drawn with
   stacked bars to visualise base composition around the motif.

4. **Integration with histone modifications (optional).**  After
   downloading histone peak BED files (e.g. from ReMap), one can intersect
   motif positions with peaks to identify motifs overlapping specific
   histone marks.  The following BEDTools command illustrates the
   workflow:

   ```bash
   bedtools intersect -a GAGG_motif_positions.bed \ 
                     -b H3K4me3_peaks.bed -wa -wb > GAGG_H3K4me3_intersect.bed
   ```

   This produces a table of motif sites that fall within H3K4me3 peaks.
   Statistics (e.g. fold enrichment) can be computed by comparing the
   number of overlapping motifs with a random expectation.

5. **Comparing assemblies (optional).**  To compare motif distributions
   across assemblies (e.g. TAIR10 vs. Col‑PEK), run the scanning script
   separately for each FASTA and then intersect or subtract the resulting
   BED files.  Differences in motif counts may highlight sequence gains or
   losses.  Whole‑genome alignment can be performed with minimap2
   (`minimap2 -ax asm5 assembly1.fa assembly2.fa > alignment.paf`) to map
   motif coordinates between assemblies.  In this environment minimap2
   could not be installed, so these commands are provided as guidance.

## Results

### Motif distribution in the Col‑CEN assembly

Scanning the Col‑CEN v1.2 assembly detected **642 988** occurrences of the
GAGG motif (and its reverse complement) across seven sequences.  The
distribution across chromosomes is summarised in Table 1.  Counts reflect
the sizes of chromosomes and the frequency of the motif; the chloroplast
(ChrC) and mitochondrion (ChrM) contain comparatively few occurrences.

**Table 1 – Number of GAGG motifs per chromosome in the Col‑CEN assembly**

| Chromosome | Motif count |
|---|---:|
| Chr1 | 156 203 |
| Chr2 | 108 340 |
| Chr3 | 128 055 |
| Chr4 | 102 490 |
| Chr5 | 144 154 |
| ChrC (chloroplast) | 740 |
| ChrM (mitochondrion) | 3 006 |

These counts show that the motif is abundant in nuclear chromosomes but rare
in organellar genomes.  This distribution may reflect differences in
nucleotide composition or epigenetic regulation between compartments.

### Sequence context around the motif

To assess base composition around the motif, a frequency matrix was
calculated from the 34‑bp windows and visualised as a stacked bar chart
(Figure 1).  The motif region (positions 15–18) shows nearly exclusive
representation of GAGG as expected, while flanking positions display a
relatively balanced composition of A, C, G and T.  Slight enrichments of
guanine immediately upstream and downstream of the motif may reflect
context preferences for 6mA modification.  The frequency matrix is saved as
a NumPy file (`motif_position_frequency.npy`) for further analyses, and the
logo plot can be regenerated or customised using the provided script.

### [Figure 1 – Nucleotide frequencies around the 6mA motif GAGG]

![Sequence logo of base frequencies around the 6mA motif]({{file:file-3rZkRPfEp4cKabUv9XUj4j}})

## Discussion

This small‑scale analysis demonstrates a pipeline for identifying 6mA
consensus motifs in a high‑quality *Arabidopsis* assembly and exploring
their local sequence context.  The motif counts reveal that GAGG sites are
common—on the order of hundreds of thousands per genome—and distributed
across all five nuclear chromosomes.  Integrating motif positions with
histone modification maps could illuminate whether 6mA preferentially occurs
within active chromatin (e.g. H3K4me3‑marked promoters) or repressive
domains (e.g. H3K9me2 heterochromatin).  The abundant histone datasets in
ReMap【404858133361378†L43-L97】 provide an opportunity for such analyses.

Comparing motif distributions between assemblies, such as TAIR10 and
Col‑PEK/Col‑CEN, may identify regions where 6mA motifs have been gained or
lost due to sequence differences.  Notably, the Col‑XJTU and Col‑CEN
assemblies add ∼14.6 Mb and 12.6 Mb of novel sequence relative to TAIR10,
including centromeric and telomeric repeats【825213852822405†L278-L297】【506180248870255†L598-L606】.  Mapping motif sites from TAIR10 onto
these assemblies using minimap2 would reveal how many motifs fall within
newly assembled regions.  Such comparative analyses could shed light on
whether 6mA is enriched in pericentromeric heterochromatin or other
structurally complex regions.

## Conclusions

N6‑methyladenine motifs are pervasive in the *Arabidopsis thaliana* genome,
with a conserved GAGG sequence context【733812400672236†L568-L575】.  Using
publicly available genome assemblies and histone datasets【404858133361378†L43-L97】,
researchers can map these motifs, compare their distribution between
assemblies and examine associations with chromatin modifications.  The
pipeline presented here demonstrates how to perform motif scanning and
generate inputs for motif discovery (MEME‑CHIP) and epigenomic integration.
Future work could expand this analysis by incorporating true histone peak
datasets and by implementing alignment‑based comparisons between TAIR10,
Col‑PEK and other assemblies to fully elucidate the regulatory role of 6mA
in plant genomes.