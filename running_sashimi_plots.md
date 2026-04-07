# Running Human MTG Sashimi Plots - Process Guide

## Summary

The script `human_MTG_sashimi.R` generates sashimi-style coverage + splice junction plots for ZDHHC8 across cell types in human MTG (middle temporal gyrus), using pre-existing merged BAM files from the Hodge et al. scRNAseq dataset.

---

## Prerequisites

### Required R Packages

| Package | Source | Purpose |
|---------|--------|---------|
| Gviz | Bioconductor | Sashimi + coverage track plotting |
| GenomicFeatures | Bioconductor | Build TxDb from GFF annotation |
| GenomicRanges | Bioconductor | Genomic interval operations |
| GenomicAlignments | Bioconductor | Read BAM junction/alignment data |
| Rsamtools | Bioconductor | BAM file I/O |
| rtracklayer | Bioconductor | Import GFF/GTF annotations |
| tidyverse | CRAN | Data wrangling (dplyr, stringr, etc.) |

### Current Cluster Situation

The cluster module `bio/R-bundle-Bioconductor/3.12-foss-2020b-R-4.0.3` has all these packages **installed** but **cannot load them** due to a missing shared library (`libssl.so.10`). This is because the cluster OS was upgraded to RHEL 9 (which ships OpenSSL 3.x) but the old R module was compiled against RHEL 7/8 (OpenSSL 1.0.x).

The newer R modules (`R/4.3.3`, `R/4.4.2`) work fine but do not have Bioconductor packages pre-installed.

---

## Setup: Create a Conda Environment (Recommended)

Since the system Bioconductor module is broken, the most reliable approach is to create a conda environment with all the required packages.

```bash
# Create environment with R and Bioconductor packages
conda create -n zdhhc8_sashimi -c conda-forge -c bioconda \
  r-base=4.3 \
  r-tidyverse \
  bioconductor-gviz \
  bioconductor-genomicfeatures \
  bioconductor-genomicranges \
  bioconductor-genomicalignments \
  bioconductor-rsamtools \
  bioconductor-rtracklayer

# Activate the environment
conda activate zdhhc8_sashimi
```

**Estimated time:** 10-15 minutes to solve and install.

### Alternative: Install into User Library with R 4.4.2

If you prefer using the system R module instead of conda:

```bash
module load R/4.4.2-gfbf-2024a

# Start R and install packages
R
```

```r
# In R:
install.packages("BiocManager")
BiocManager::install(c("Gviz", "GenomicFeatures", "GenomicRanges",
                        "GenomicAlignments", "Rsamtools", "rtracklayer"))
install.packages("tidyverse")
```

This installs to `~/R/x86_64-pc-linux-gnu-library/4.4/`. Takes ~20-30 minutes.

---

## Running the Script

### Option A: Interactive R Session

```bash
# If using conda:
conda activate zdhhc8_sashimi

# If using module:
module load R/4.4.2-gfbf-2024a

# Run the script
cd /external/rprshnas01/kcni/stripathy/isoform_project
Rscript human_MTG_sashimi.R
```

### Option B: Submit as a Cluster Job

```bash
cat > /external/rprshnas01/kcni/stripathy/isoform_project/run_sashimi.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=zdhhc8_sashimi
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=sashimi_%j.out

# Use whichever R setup you chose:
# conda activate zdhhc8_sashimi
# OR
# module load R/4.4.2-gfbf-2024a

cd /external/rprshnas01/kcni/stripathy/isoform_project
Rscript human_MTG_sashimi.R
EOF

sbatch run_sashimi.sh
```

### Option C: RStudio Server

```bash
module load RStudio-Server/2023.12.1+402-gfbf-2023b-Java-11-R-4.3.3
```

Open `human_MTG_sashimi.R` in RStudio and source it. Note: you'll still need Bioconductor packages installed (see setup above).

---

## What the Script Does

1. **Builds a TxDb** from the human RefSeq GRCh38.p2 GFF annotation (the same one used for the original STAR alignment). This is the slowest step (~3-5 minutes).

2. **Defines `makefigure_human()`** — adapted from Mel's `makefigure()` in `sample_viz.Rmd`. Key changes from the mouse version:
   - Chromosome: `NC_000022.11` (human chr22) instead of `NC_000082.7` (mouse chr16)
   - Junction coordinates: start=20143757, end=20145228 (short/OV) or 20147020 (long/EV) instead of the mouse coordinates
   - Gene range: 20131841-20148007

3. **Generates 4 TIFF figures** saved to `human_MTG_figures/`:

   | File | Contents |
   |------|----------|
   | `MTG_exc_inh_astro.tiff` | Excitatory vs Inhibitory vs Astrocyte (3 panels) |
   | `MTG_exc_subclasses_1.tiff` | L5-6 IT Car3, L5-6 NP, L6 CT, L6b (4 panels) |
   | `MTG_inh_subclasses.tiff` | PVALB, SST, LAMP5, VIP (4 panels) |
   | `MTG_exc_vs_inh.tiff` | Excitatory vs Inhibitory (2 panels) |

Each figure shows:
- **Top track:** Genomic coordinate axis
- **Middle tracks:** Read coverage (filled area) + sashimi arcs (junction-spanning reads) per cell type BAM
- **Bottom track:** ZDHHC8 gene model with transcript variant annotations (NM_001185024.1 and NM_013373.3)

---

## Input Data

All input BAMs are pre-built merged BAMs at:
```
/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/Isoform_project/Hodge_aligned/MTG/merged_bams/
```

These were created by Mel by:
1. Downloading Hodge et al. scRNAseq FASTQs from NEMO (Allen Institute)
2. Aligning with STAR to RefSeq GRCh38.p2
3. Randomly sampling 100 cells per subclass
4. Merging per-cell BAMs into subclass-level BAMs, filtered to the ZDHHC8 locus

See `data_inventory.md` for the full data map.

---

## Extending to Other Brain Regions

The same script can be adapted for **V1** and **S1** by changing the BAM directory path:

```r
# For V1 (visual cortex):
bam_dir <- ".../Hodge_aligned/V1/merged_bams"

# For S1 (somatosensory cortex):
bam_dir <- ".../Hodge_aligned/S1/merged_bams"
```

V1 has the same cell types as MTG. S1 has fewer subclass-level BAMs (mainly excitatory and inhibitory pooled).

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| `libssl.so.10: cannot open shared object` | Don't use the old Bioconductor module. Use conda or R 4.4 with user-installed packages. |
| `makeTxDbFromGFF` takes very long | The GFF is large (~1.5 GB). This is normal; it takes 3-5 min. Only runs once per session. |
| Empty sashimi arcs | The merged BAMs are small (filtered to ZDHHC8 locus only). Some rare cell types may have very few junction-spanning reads. The script sets a minimum threshold of 5. |
| `Error in readGAlignmentPairs` | Make sure the `.bam.bai` index file exists alongside each BAM. They should already be present. |
