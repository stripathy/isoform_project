# ZDHHC8 Isoform Project - Data Inventory

## Overview

This document catalogs the data and code locations for the ZDHHC8 isoform splicing project, a collaboration between the Tripathy Lab (CAMH/U of T) and the Bamji Lab (UBC). The project examines cell-type-specific alternative splicing of ZDHHC8 in mouse and human brain using scRNAseq data.

Former RA: Melanie Davies (GitHub: daviemel)

---

## Code

### GitHub Repository
- **URL:** https://github.com/daviemel/isoform_project
- **Local clone:** `/external/rprshnas01/kcni/stripathy/isoform_project/`
- **Key files:**
  - `sample_viz.Rmd` — Main notebook: FASTQ sampling, STAR alignment, BAM merging, sashimi plot generation (mouse data)
  - `isoform-viz-functions.R` — Refactored/reusable versions of sashimi plot and exon coordinate functions
  - `zdhhc8-isoform-quant-functions.R` — BAM indexing, merging, junction-based and read-based splicing ratio quantification (supports both mouse and human coordinates)
  - `STAR_align.sh` — Shell script to generate per-sample STAR alignment commands from FASTQ files
  - `gtex_analysis/gtex_zdhhc8_junction_analysis.R` — GTEx bulk tissue junction ratio analysis
  - `gtex_analysis/zdhhc8_junction_data_subsetting.R` — Subsets the 14 GB GTEx junction file to just ZDHHC8

### Local GTEx Analysis (Shreejoy)
- **Location:** `/external/rprshnas01/kcni/stripathy/zdhhc8_splicing/`
- `scratch.R` — GTEx ZDHHC8 junction ratio boxplots by tissue (broad, fine, brain-only)
- `zdhhc8_junction_subsetting.R` — Same GTEx junction subsetting as in GitHub repo
- `data/zdhhc8_junctions.csv` — Pre-extracted ZDHHC8 junction counts from GTEx v8

---

## Human Data (Hodge et al. scRNAseq)

### Source
Human MTG, V1, and S1 scRNAseq from the Allen Institute (Hodge et al.), downloaded via NEMO:
`/external/rprshnas01/netdata_kcni/stlab/NEMO/lein-human-cortex/bundled/transcriptome/scell/SSv4/human_gru/raw/`

### Reference Genome
- **Assembly:** RefSeq GRCh38.p2
- **GFF:** `/external/rprshnas01/netdata_kcni/stlab/Genomic_references/AIBS/Human/Refseq_GRCh38.p2/Raw/GCF_000001405.28_GRCh38.p2_genomic.gff`
- **STAR index:** `/external/rprshnas01/netdata_kcni/stlab/Genomic_references/AIBS/Human/Refseq_GRCh38.p2/USE_THIS_genomeDir_gff/`
- **ZDHHC8 locus:** NC_000022.11:20131841-20148007
- **Transcripts:** NM_001185024.1 (variant 1, long/EV), NM_013373.3 (variant 2, short/OV)

### STAR-Aligned Data
**Base path:** `/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/Isoform_project/Hodge_aligned/`

#### MTG (Middle Temporal Gyrus)
- **Cell type directories:** Astrocyte, PVALB, SST, VIP, LAMP5, L5-6_IT_Car3, L5-6_NP, L6_CT, L6b (each has `_outputs/STAR_results/`)
- **Merged BAMs (ready to use):** `MTG/merged_bams/`
  - Per-subclass: `Astrocyte.bam`, `PVALB.bam`, `SST.bam`, `VIP.bam`, `LAMP5.bam`, `L5-6_IT_Car3.bam`, `L5-6_NP.bam`, `L6_CT.bam`, `L6b.bam`
  - Per-class: `excitatory.bam`, `inhibitory.bam`
  - All have corresponding `.bam.bai` index files

#### V1 (Visual Cortex)
- **Same cell types as MTG**
- **Merged BAMs:** `V1/merged_bams/` — same structure as MTG (per-subclass + excitatory/inhibitory)

#### S1 (Somatosensory Cortex)
- **Cell type directories:** IT, L5-6_IT_Car3, L5-6_NP, L6b, SST, VIP, LAMP5
- **Merged BAMs:** `S1/merged_bams/` — `excitatory.bam` (split into `excitatory_1.bam` through `excitatory_4.bam`), `inhibitory.bam`
- Note: S1 has fewer subclass-level merged BAMs than MTG/V1

---

## Mouse Data (AIBS scRNAseq 2019 / Yao et al.)

### Source
- **Raw FASTQs:** `/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_scRNAseq_2019/mouse/`
- **Metadata:** `/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_scRNAseq_2019/mouse/2020_10_24/metadata.csv`

### Reference Genome
- **Assembly:** RefSeq GRCm39
- **GFF:** `/external/rprshnas01/netdata_kcni/stlab/Genomic_references/Refseq/Mouse/Refseq_GRCm39/Raw/GCF_000001635.27_GRCm39_genomic.gff`
- **STAR index:** `/external/rprshnas01/netdata_kcni/stlab/Genomic_references/Refseq/Mouse/Refseq_GRCm39/USE_THIS_genomeDir_gff`
- **Zdhhc8 locus:** NC_000082.7:18038612-18056471

### Processed Data
- **Salmon quantification:** `/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/Isoform_project/Yao_salmon_quant/`
  - Contains ZDHHC8 isoform counts for ACA and VISp (`visp_zdhhc8_isoform_counts.csv`, `aca_zdhhc8_isoform_counts.csv`)
  - Also has DHODH and KCN isoform counts
- **Random 100 cell sampling (Yao):** `Yao_Random100/` — SSp and HIP directories, appear mostly empty
- **Note:** Mouse merged BAMs from the Yao dataset are not present on disk; they were generated on-the-fly per the notebook in the GitHub repo

---

## Human Bulk Tissue Data (GTEx v8)

### Source
- **GTEx junction file:** `/external/rprshnas01/netdata_kcni/stlab/Public/GTEx/GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.gct` (14 GB)
- **Pre-extracted ZDHHC8 junctions:** `/external/rprshnas01/netdata_kcni/stlab/Public/GTEx/zdhhc8_junctions.csv`
- **GTEx sample metadata:** Downloaded at runtime from `https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt`

### Key Junction Identifiers
- Excitatory-specific (long isoform): `chr22_20143757_20147020`
- Inhibitory-specific (short isoform): `chr22_20143757_20145228`
- Splicing ratio: `exc_spec / (exc_spec + int_spec)`

---

## What's Missing or Incomplete

- **Interspecies comparison directory** (`Hodge_aligned/Interspecies_comparison/`) — has `human/` and `mouse/` subdirectories but they appear empty
- **Yao Random100** — SSp and HIP directories exist but contain no data files
- **Mouse merged BAMs** — not saved to disk; would need to re-run the BAM merging pipeline from `sample_viz.Rmd`
- **No local clone of the GitHub repo existed** prior to this investigation (now cloned to `~/isoform_project/`)
