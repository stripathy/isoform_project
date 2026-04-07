# Mouse Zdhhc8 Sashimi Plot & Junction Quantification Pipeline

## Goal

Reproduce the human MTG ZDHHC8 sashimi plot and junction quantification analysis for **mouse brain data** (Yao et al. AIBS scRNAseq 2019). Specifically:

1. Align Smart-seq2 BAM/FASTQ files from mouse single-cell data to the **Zdhhc8 locus** on chr16
2. Merge per-cell BAMs into combined **excitatory** and **inhibitory** subclass-level BAMs
3. Generate **sashimi plots** showing read coverage and splice junction arcs per cell type
4. **Quantify junction usage** supporting isoform 1 (long/EV) vs isoform 2 (short/OV) per cell type population

---

## Pipeline Overview

```
Raw FASTQs (per cell)
    │
    ▼
[STAR alignment] ──► Per-cell coordinate-sorted BAMs
    │
    ▼
[makeSubclassBAMs()] ──► Subclass-level merged BAMs (filtered to Zdhhc8 locus)
    │
    ├──► [Sashimi plots via Gviz]        ──► SVG/TIFF figures
    │
    └──► [getJunctions() / getReads()]    ──► Isoform ratio quantification (CSV)
```

---

## Step 1: STAR Alignment

### Recommended: Locus-Only Alignment (Fast)

For feasibility testing (and likely for the final analysis), you can build a **Zdhhc8-locus-only STAR index** instead of using the full genome. This reduces the index from ~30 GB to a few MB and makes alignment near-instant per cell.

#### 1a. Extract the Zdhhc8 Region from the Genome FASTA

Extract chromosome 16 (or a padded region around Zdhhc8) from the full genome FASTA. Adding ±50 kb flanking ensures splice junction reads near the locus edges are captured.

```bash
# First, index the genome FASTA if not already done
samtools faidx /path/to/GCF_000001635.27_GRCm39_genomic.fna

# Option A: Extract all of chr16 (simple, still much smaller than full genome)
samtools faidx /path/to/GCF_000001635.27_GRCm39_genomic.fna NC_000082.7 > zdhhc8_chr16.fa

# Option B: Extract just the Zdhhc8 locus ± 50 kb flanking (smallest possible)
samtools faidx /path/to/GCF_000001635.27_GRCm39_genomic.fna \
  NC_000082.7:17988612-18106471 > zdhhc8_region.fa
```

#### 1b. Subset the GFF Annotation to the Region

Only include annotation features that overlap the extracted region. This avoids warnings about missing chromosomes and speeds up TxDb construction later.

```bash
# For chr16-only:
grep -E "^#|^NC_000082\.7" /path/to/GCF_000001635.27_GRCm39_genomic.gff > zdhhc8_chr16.gff

# For the locus-only region, further filter by coordinate range:
awk -F'\t' 'BEGIN{OFS="\t"} /^#/{print; next} $1=="NC_000082.7" && $4 <= 18106471 && $5 >= 17988612 {print}' \
  /path/to/GCF_000001635.27_GRCm39_genomic.gff > zdhhc8_region.gff
```

#### 1c. Build a Small STAR Index

```bash
mkdir zdhhc8_star_index

STAR --runMode genomeGenerate \
     --genomeDir zdhhc8_star_index \
     --genomeFastaFiles zdhhc8_chr16.fa \
     --sjdbGTFfile zdhhc8_chr16.gff \
     --sjdbGTFtagExonParentTranscript Parent \
     --genomeSAindexNbases 10 \
     --runThreadN 4
```

**Notes:**
- `--genomeSAindexNbases 10` is required for small genomes/regions (default 14 is too large). STAR will suggest the optimal value in its error message if this needs adjustment.
- `--sjdbGTFtagExonParentTranscript Parent` is needed for GFF3 files (vs GTF which uses `transcript_id`).
- This should complete in under a minute.

#### 1d. Align Using the Locus-Only Index

Use the same `STAR_align.sh` script, just point `genomePath` to the new small index:

```bash
bash STAR_align.sh \
  /path/to/fastq_dir \
  /path/to/output_dir \
  /path/to/zdhhc8_star_index \
  /path/to/zdhhc8_chr16.gff \
  /path/to/filesofInterest.csv
```

Each cell should align in seconds rather than minutes.

#### Caveat

With a locus-only index, reads that truly belong elsewhere in the genome have no competing alignment target, so a small number may falsely map to the Zdhhc8 region. In practice this is negligible for junction-level analysis because:
- False mappings are unlikely to span the specific isoform-defining splice junctions
- The sashimi plot coverage signal is dominated by true mappings
- Junction quantification filters for exact start/end positions

For a paranoia check, you could align a subset of cells to the full genome and compare junction counts — they should be nearly identical.

---

### Alternative: Full-Genome Alignment (Standard)

`STAR_align.sh` takes paired-end Smart-seq2 FASTQs and aligns them to the full reference genome. More accurate but much slower and requires the full ~30 GB STAR index. This is what the original pipeline used.

### Usage

```bash
bash STAR_align.sh <datadir> <outputdir> <genomePath> <annotationPath> <filesofInterest.csv>
```

| Argument | Description |
|----------|-------------|
| `datadir` | Directory containing `.fastq.tar` archives (one per cell) |
| `outputdir` | Where to write outputs |
| `genomePath` | Path to pre-built STAR genome index directory |
| `annotationPath` | Path to genome annotation (GFF/GTF) |
| `filesofInterest.csv` | CSV with sample names (one per line, no extension) |

### STAR Parameters

```
--runThreadN 12
--sjdbGTFfile <annotation>
--quantMode GeneCounts
--outSAMtype BAM SortedByCoordinate
--readFilesCommand zcat
```

### Expected Input File Structure

Each sample in the CSV should have a corresponding `{samplename}.fastq.tar` in the data directory, containing:
- `{samplename}_R1.fastq.gz`
- `{samplename}_R2.fastq.gz`

### Output Structure

```
<outputdir>/
├── STAR_scripts/           # One .sh script per sample (for parallel submission)
│   ├── STARParaScript0.sh
│   ├── STARParaScript1.sh
│   └── ...
├── STARParaCom.txt         # List of all scripts (for batch submission)
└── STAR_results/
    ├── coord_bams/         # Coordinate-sorted BAMs (one per cell)
    ├── pcrless_bams/       # (unused in current pipeline)
    └── star_logs/
        ├── Init_QCMetrics/ # STAR Log.final.out per sample
        └── Init_Runmsgs/   # STAR Log.out per sample
```

### Cluster Submission

The script generates individual shell scripts for each sample. Submit them in parallel:

```bash
# Option A: GNU parallel
cat STARParaCom.txt | parallel -j 10 bash {}

# Option B: SLURM array job
# Wrap each script line in an sbatch submission
while read script; do
    sbatch --mem=32G --cpus-per-task=12 --time=01:00:00 "$script"
done < STARParaCom.txt
```

### Required Data

| Item | Path on Cluster | Size |
|------|----------------|------|
| Mouse STAR genome index | `/external/rprshnas01/netdata_kcni/stlab/Genomic_references/Refseq/Mouse/Refseq_GRCm39/USE_THIS_genomeDir_gff` | ~30 GB |
| Mouse GFF annotation | `/external/rprshnas01/netdata_kcni/stlab/Genomic_references/Refseq/Mouse/Refseq_GRCm39/Raw/GCF_000001635.27_GRCm39_genomic.gff` | ~1.5 GB |
| Raw FASTQs (Yao et al.) | `/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_scRNAseq_2019/mouse/` | Variable |
| Cell metadata | `/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_scRNAseq_2019/mouse/2020_10_24/metadata.csv` | — |

---

## Step 2: Merge Per-Cell BAMs into Subclass BAMs

### What It Does

`makeSubclassBAMs()` in `zdhhc8-isoform-quant-functions.R`:

1. For each subclass, finds all per-cell BAMs in `{bamdir}/{subclass}_outputs/STAR_results/coord_bams/`
2. **Randomly samples 100 cells** per subclass (seed = 12345)
3. Merges the 100 BAMs into one, **filtered to the Zdhhc8 gene locus only**
4. Indexes the output BAM

### Usage (in R)

```r
source("zdhhc8-isoform-quant-functions.R")

# Build TxDb from the mouse GFF
txdb <- makeTxDbFromGFF("/path/to/GCF_000001635.27_GRCm39_genomic.gff")

# Define subclasses (should match folder names in bamdir minus "_outputs")
exc_subclasses <- c("L5_IT", "L6_IT", "L6_CT", "L6b", "NP")  # adjust to actual folder names
inh_subclasses <- c("Pvalb", "Sst", "Lamp5", "Vip")            # adjust to actual folder names

# Create merged BAMs
makeSubclassBAMs(
  bamdir = "/path/to/aligned/VISp",
  outputdir = "/path/to/output/merged_bams",
  subclasses = exc_subclasses,
  txdb = txdb,
  gene = "Zdhhc8",
  index = TRUE  # set TRUE if per-cell BAMs are not yet indexed
)
```

### Creating Class-Level BAMs (Excitatory / Inhibitory)

After subclass BAMs are created, merge them into class-level BAMs using samtools:

```bash
# Excitatory = merge all excitatory subclass BAMs
samtools merge excitatory.bam L5_IT.bam L6_IT.bam L6_CT.bam L6b.bam NP.bam
samtools index excitatory.bam

# Inhibitory = merge all inhibitory subclass BAMs
samtools merge inhibitory.bam Pvalb.bam Sst.bam Lamp5.bam Vip.bam
samtools index inhibitory.bam
```

### Key Detail: 100-Cell Sampling

The sampling happens at the BAM merge step, not during alignment. Each subclass contributes exactly 100 randomly sampled cells, so:
- Excitatory BAM ≈ 400-500 cells (4-5 subclasses × 100)
- Inhibitory BAM ≈ 400 cells (4 subclasses × 100)

---

## Step 3: Junction Quantification

### What It Does

`getJunctions()` counts reads spanning the **isoform-defining splice junction** to compute isoform ratios per cell type.

### How the Two Isoforms Are Distinguished

Both isoforms share a common splice junction **start** position. They differ in where that junction **ends**, reflecting exon inclusion vs. skipping:

**Human coordinates (chr22, NC_000022.11):**
- Junction start: **20,143,757**
- Long isoform (EV, NM_001185024.1): junction end at **20,147,020**
- Short isoform (OV, NM_013373.3): junction end at **20,145,228**

**Mouse coordinates (chr16, NC_000082.7):**
- Zdhhc8 locus: **18,038,612 – 18,056,471**
- **Mouse-specific junction coordinates need to be identified from the GFF.** The human coordinates above are hardcoded in `getJunctions()` and `getReads()` — these must be updated with the orthologous mouse positions.

### Identifying Mouse Junction Coordinates

Run this on the cluster to extract mouse Zdhhc8 transcript structure from the GFF:

```r
library(GenomicFeatures)
library(tidyverse)

txdb <- makeTxDbFromGFF("/path/to/GCF_000001635.27_GRCm39_genomic.gff")

# Find Zdhhc8 gene boundaries
genes_txdb <- genes(txdb, single.strand.genes.only = FALSE)
gene_info <- data.frame(genes_txdb) %>%
  filter(group_name == "Zdhhc8" & str_starts(seqnames, "NC"))
print(gene_info)

# List all transcripts in the locus
chrID <- as.character(gene_info[,3])
coords <- c(gene_info[,4], gene_info[,5])

tx <- data.frame(transcripts(txdb)) %>%
  filter(seqnames == chrID) %>%
  filter(start >= coords[1] & end <= coords[2]) %>%
  filter(str_starts(tx_name, "NM") | str_starts(tx_name, "XM"))
print(tx)

# Get exon boundaries for each transcript
exons <- exonsBy(txdb, by = "tx", use.names = TRUE)[tx$tx_name]
print(as.data.frame(exons))
```

Look for the two NM_ transcripts and identify the exon boundary that differs between them. The junction start = end of the shared upstream exon + 1, and the two junction ends = start of the next exon - 1 for each isoform.

### Usage (After Identifying Mouse Coordinates)

You will need to update the junction coordinates in `getJunctions()` and `getReads()` for mouse, or parameterize them. Then:

```r
source("zdhhc8-isoform-quant-functions.R")

# List merged BAMs and their subclass labels
bam_files <- c(
  "/path/to/merged_bams/Pvalb.bam",
  "/path/to/merged_bams/Sst.bam",
  "/path/to/merged_bams/Lamp5.bam",
  "/path/to/merged_bams/Vip.bam",
  "/path/to/merged_bams/L5_IT.bam",
  "/path/to/merged_bams/L6_CT.bam"
)
subclass_names <- c("Pvalb", "Sst", "Lamp5", "Vip", "L5_IT", "L6_CT")

# Junction-based quantification
junction_results <- getJunctions(bam_files, subclass_names)
print(junction_results)
write.csv(junction_results, "mouse_zdhhc8_junction_ratios.csv", row.names = FALSE)

# Read-based quantification (alternative validation)
read_results <- getReads(bam_files, subclass_names)
print(read_results)
write.csv(read_results, "mouse_zdhhc8_read_ratios.csv", row.names = FALSE)
```

### Output

Both functions return a dataframe with columns:

| subclass | EV (%) | OV (%) |
|----------|--------|--------|
| Pvalb | 35.2 | 64.8 |
| Sst | 28.1 | 71.9 |
| L5_IT | 72.4 | 27.6 |
| ... | ... | ... |

Where EV = excitatory variant (long isoform) and OV = other variant (short isoform).

---

## Step 4: Sashimi Plot Visualization

### What It Does

Uses Gviz to create publication-quality sashimi plots showing:
- **Coverage tracks** (filled area plots) per cell type BAM
- **Sashimi arcs** connecting splice junctions, with arc thickness proportional to read count
- **Gene model track** showing exon/intron structure for all NM_ transcripts

### Setup: R Environment on the Cluster

The cluster's old Bioconductor module is broken (missing `libssl.so.10`). Use conda:

```bash
conda create -n zdhhc8_sashimi -c conda-forge -c bioconda \
  r-base=4.3 \
  r-tidyverse \
  r-svglite \
  bioconductor-gviz \
  bioconductor-genomicfeatures \
  bioconductor-genomicranges \
  bioconductor-genomicalignments \
  bioconductor-rsamtools \
  bioconductor-rtracklayer

conda activate zdhhc8_sashimi
```

Alternatively, install into a user library with R 4.4.2:

```bash
module load R/4.4.2-gfbf-2024a
Rscript -e 'install.packages("BiocManager"); BiocManager::install(c("Gviz", "GenomicFeatures", "GenomicRanges", "GenomicAlignments", "Rsamtools", "rtracklayer")); install.packages(c("tidyverse", "svglite"))'
```

### Adapting the Script for Mouse

The existing `human_MTG_sashimi.R` can be adapted. Key changes needed:

| Parameter | Human | Mouse |
|-----------|-------|-------|
| Chromosome | `NC_000022.11` | `NC_000082.7` |
| Gene locus | `20131841-20148007` | `18038612-18056471` |
| GFF path | `GCF_000001405.28_GRCh38.p2_genomic.gff` | `GCF_000001635.27_GRCm39_genomic.gff` |
| Junction start | `20143757` | **TBD from GFF** |
| EV junction end | `20147020` | **TBD from GFF** |
| OV junction end | `20145228` | **TBD from GFF** |

### Generating Figures

```r
# Example: adapted makefigure for mouse (pseudocode)
svglite("mouse_VISp_exc_vs_inh.svg", width = 10, height = 8)
makefigure_mouse(
  txdb,
  bam1 = "merged_bams/excitatory.bam",
  bam2 = "merged_bams/inhibitory.bam",
  chrID = "NC_000082.7",
  coords = c(18038612, 18056471),
  bam1title = "Excitatory",
  bam2title = "Inhibitory"
)
dev.off()
```

Use `svglite()` instead of `svg()` — the base R `svg()` device requires Cairo/X11, which may not be available on the cluster. Install `svglite` via conda or `install.packages("svglite")`.

---

## Feasibility Testing Checklist

### Minimal Test: Can We Build the TxDb and Find Mouse Zdhhc8?

**Requires:** Mouse GFF only (no BAMs needed)

```r
txdb <- makeTxDbFromGFF("/path/to/GCF_000001635.27_GRCm39_genomic.gff")
genes_txdb <- genes(txdb, single.strand.genes.only = FALSE)
gene_info <- data.frame(genes_txdb) %>% filter(group_name == "Zdhhc8")
print(gene_info)  # should show NC_000082.7:18038612-18056471
```

**Purpose:** Confirms the GFF is readable and the gene lookup works. Also use this to extract the mouse-specific junction coordinates (see Step 3 above).

### Small-Scale Test: Align a Few Cells and Make a Sashimi Plot

1. Pick ~10 cells from one excitatory and one inhibitory subclass from the metadata
2. Run STAR alignment for just those 20 cells
3. Merge into two small BAMs (exc/inh) filtered to Zdhhc8
4. Generate a sashimi plot

**Requires:** STAR genome index + GFF + 20 FASTQ archives + metadata

### Full Run: All Cell Types, Junction Quantification

1. Identify all cells per subclass from the metadata CSV
2. Align all cells with STAR (batch job)
3. Run `makeSubclassBAMs()` for all subclasses (100 cells each)
4. Merge into class-level BAMs
5. Run `getJunctions()` and `getReads()`
6. Generate full sashimi figure set

**Requires:** Everything above at full scale

---

## Data Locations on Cluster

| Item | Path |
|------|------|
| **Mouse FASTQs** | `/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_scRNAseq_2019/mouse/` |
| **Cell metadata** | `/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_scRNAseq_2019/mouse/2020_10_24/metadata.csv` |
| **Mouse GFF** | `/external/rprshnas01/netdata_kcni/stlab/Genomic_references/Refseq/Mouse/Refseq_GRCm39/Raw/GCF_000001635.27_GRCm39_genomic.gff` |
| **STAR genome index** | `/external/rprshnas01/netdata_kcni/stlab/Genomic_references/Refseq/Mouse/Refseq_GRCm39/USE_THIS_genomeDir_gff` |
| **Code repo** | `/external/rprshnas01/kcni/stripathy/isoform_project/` (or clone from GitHub) |
| **Existing Salmon quant** | `/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/Isoform_project/Yao_salmon_quant/` |

### Brain Regions Available in Yao et al.

- ACA (Anterior Cingulate Area)
- VISp (Primary Visual Cortex)
- SSp (Primary Somatosensory Cortex)
- HIP (Hippocampus)

---

## Key Code Files

| File | Purpose |
|------|---------|
| `STAR_align.sh` | Generates per-sample STAR alignment scripts from FASTQs |
| `zdhhc8-isoform-quant-functions.R` | `indexBAMs()`, `makeSubclassBAMs()`, `getJunctions()`, `getReads()` |
| `isoform-viz-functions.R` | `exonCoords()`, `junctionsGRanges()`, `coverageJunctionPlot()` — generic sashimi plotting |
| `human_MTG_sashimi.R` | Human-specific sashimi script (template for mouse adaptation) |
| `sample_viz.Rmd` | Original notebook by Mel covering the full mouse pipeline end-to-end |

---

## Important Notes

- **Junction coordinates are hardcoded** in `getJunctions()` and `getReads()` for human (start=20143757, EV end=20147020, OV end=20145228). These **must be updated** with mouse-orthologous positions before running quantification on mouse data.
- **`makeSubclassBAMs()` samples exactly 100 cells** per subclass with a fixed seed (12345). If a subclass has fewer than 100 cells, the function will error — add a check or adjust the sample size.
- **Mouse merged BAMs were not saved** from previous runs. They need to be regenerated from scratch.
- The `sample_viz.Rmd` notebook contains the original mouse pipeline and may be useful as a reference, but the R functions have since been refactored into the standalone `.R` files listed above.
