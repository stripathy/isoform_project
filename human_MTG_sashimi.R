##############################################################################
# Human MTG Sashimi Plots for ZDHHC8
# Adapted from sample_viz.Rmd (daviemel/isoform_project)
#
# Uses merged BAMs from Hodge et al. human MTG scRNAseq data,
# aligned to RefSeq GRCh38.p2
##############################################################################

# Load Packages
suppressPackageStartupMessages({
  library(Gviz)
  library(rtracklayer)
  library(Rsamtools)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(GenomicAlignments)
  library(tidyverse)
  library(svglite)
})

# =============================================================================
# 1. Reference Genome & TxDb
# =============================================================================

# Human RefSeq GRCh38.p2 GFF (same annotation used for STAR alignment)
human_gff <- file.path(getwd(), "data", "GCF_000001405.28_GRCh38.p2_genomic.gff")

message("Building TxDb from human GFF (this may take several minutes)...")
humanTxDb <- makeTxDbFromGFF(human_gff)
message("TxDb ready.")

# =============================================================================
# 2. Human ZDHHC8 Coordinates
# =============================================================================

# ZDHHC8 on chr22 (RefSeq accession NC_000022.11)
# Gene range: 20131841 - 20148007
# Transcripts:
#   NM_001185024.1 = transcript variant 1 (long / "EV" = excitatory variant)
#   NM_013373.3    = transcript variant 2 (short / "OV" = other variant)
#
# Key splice junctions (from zdhhc8-isoform-quant-functions.R):
#   Both start at 20143757
#   OV junction end: 20145228
#   EV junction end: 20147020

chrID_human <- "NC_000022.11"
coords_human <- c(20131841, 20148007)

# Junction coordinates for sashimi score filtering
jxn_start <- 20143757
jxn_end_ov <- 20145228   # short isoform junction
jxn_end_ev <- 20147020   # long isoform junction

# =============================================================================
# 3. Merged BAM Paths (Human MTG - Hodge et al.)
# =============================================================================

bam_dir <- file.path(getwd(), "data", "MTG", "merged_bams")

# Class-level BAMs
excitatory <- file.path(bam_dir, "excitatory.bam")
inhibitory <- file.path(bam_dir, "inhibitory.bam")

# Subclass-level BAMs
astrocyte    <- file.path(bam_dir, "Astrocyte.bam")
l56_it_car3  <- file.path(bam_dir, "L5-6_IT_Car3.bam")
l56_np       <- file.path(bam_dir, "L5-6_NP.bam")
l6_ct        <- file.path(bam_dir, "L6_CT.bam")
l6b          <- file.path(bam_dir, "L6b.bam")
lamp5        <- file.path(bam_dir, "LAMP5.bam")
pvalb        <- file.path(bam_dir, "PVALB.bam")
sst          <- file.path(bam_dir, "SST.bam")
vip          <- file.path(bam_dir, "VIP.bam")

# =============================================================================
# 4. Sashimi Plot Function (adapted for human ZDHHC8)
# =============================================================================

makefigure_human <- function(txdb, bam1, bam2, bam3 = NULL, bam4 = NULL,
                             chrID, coords,
                             bam1title = deparse(substitute(bam1)),
                             bam2title = deparse(substitute(bam2)),
                             bam3title = deparse(substitute(bam3)),
                             bam4title = deparse(substitute(bam4))) {

  # Extract list of transcripts in the region
  gene_tx <- data.frame(transcripts(txdb)) %>%
    filter(seqnames == chrID) %>%
    filter((start >= coords[1]) & (end <= coords[2])) %>%
    dplyr::select(tx_name) %>%
    dplyr::filter(str_starts(tx_name, "NM")) %>%
    pull(tx_name)

  # Create genome axis track
  genomeAxis <- GenomeAxisTrack(
    name = "MyAxis",
    col = "lightsteelblue4",
    fontcolor = "lightsteelblue4",
    add35 = TRUE
  )

  bams <- c(bam1, bam2, bam3, bam4)
  bamtitles <- c(bam1title, bam2title, bam3title, bam4title)

  # Prepare plot for each BAM
  i <- 1
  bamplots <- vector()

  for (bam in bams) {
    if (!is.null(bam)) {
      # Read alignment pairs and get junction info
      galign <- readGAlignmentPairs(file = bam, index = bam, strandMode = 1)
      junctions <- summarizeJunctions(galign)
      junction_df <- data.frame(
        c(as.data.frame(junctions@ranges), as.data.frame(junctions$score))
      )

      # Filter for the two ZDHHC8-specific junctions (human coordinates)
      junction_df <- junction_df %>%
        filter(start == jxn_start) %>%
        filter(end == jxn_end_ov | end == jxn_end_ev)

      # Compute min junction threshold so both arcs are visible
      junction_min <- min(junction_df$junctions.score, na.rm = TRUE)

      if (is.infinite(junction_min) | junction_min < 5 | is.na(junction_min)) {
        junction_min <- 5
      }

      # Create coverage + sashimi track
      bamplot <- AlignmentsTrack(
        bam,
        name = bamtitles[i],
        cex = 2,
        background.title = "white",
        sashimiScore = junction_min,
        col.axis = "lightsteelblue4",
        col.title = "lightsteelblue4"
      )
      bamplots <- append(bamplots, bamplot)
      i <- i + 1
    }
  }

  # Create gene model track
  gr <- exonsBy(txdb, by = "tx", use.names = TRUE)[gene_tx]
  gr <- unlist(gr)
  elementMetadata(gr)$transcript <- names(gr)
  gene_models <- Gviz::GeneRegionTrack(
    gr,
    showId = TRUE,
    options(ucscChromosomeNames = FALSE),
    transcriptAnnotation = "transcript",
    name = "Gene Model",
    background.title = "white",
    col.axis = "lightsteelblue4",
    col.title = "lightsteelblue4",
    fill = "darkgrey",
    fontcolor.group = "lightsteelblue4",
    col.line = "lightsteelblue4"
  )

  # Assemble tracks
  tracks <- c(genomeAxis)
  for (plot in bamplots) {
    tracks <- append(tracks, plot)
  }
  tracks <- append(tracks, gene_models)

  # Track sizes
  tracksizes <- c(1)
  for (bam in bams) {
    if (!is.null(bam)) {
      tracksizes <- append(tracksizes, 3)
    }
  }
  tracksizes <- append(tracksizes, 2)

  # Render figure
  options(ucscChromosomeNames = FALSE)
  fig <- plotTracks(
    trackList = tracks,
    showId = TRUE,
    transcriptAnnotation = "transcript",
    chromosome = chrID,
    sizes = tracksizes,
    from = coords[1],
    to = coords[2],
    extend.left = 3500,
    fill = "lightsteelblue",
    col.sashimi = "lightsteelblue4",
    type = c('coverage', 'sashimi')
  )

  return(fig)
}

# =============================================================================
# 5. Generate Sashimi Plots
# =============================================================================

output_dir <- file.path(getwd(), "human_MTG_figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Plot 1: Excitatory vs Inhibitory vs Astrocyte ---
message("Generating: Excitatory vs Inhibitory vs Astrocyte...")
svglite(file.path(output_dir, "MTG_exc_inh_astro.svg"),
     width = 10, height = 8)
makefigure_human(
  humanTxDb,
  bam1 = excitatory, bam2 = inhibitory, bam3 = astrocyte,
  chrID = chrID_human, coords = coords_human,
  bam1title = "Excitatory", bam2title = "Inhibitory", bam3title = "Astrocyte"
)
dev.off()

# --- Plot 2: Excitatory subclasses (L5-6 IT Car3, L5-6 NP, L6 CT, L6b) ---
message("Generating: Excitatory subclasses (group 1)...")
svglite(file.path(output_dir, "MTG_exc_subclasses_1.svg"),
     width = 10, height = 10)
makefigure_human(
  humanTxDb,
  bam1 = l56_it_car3, bam2 = l56_np, bam3 = l6_ct, bam4 = l6b,
  chrID = chrID_human, coords = coords_human,
  bam1title = "L5-6 IT Car3", bam2title = "L5-6 NP",
  bam3title = "L6 CT", bam4title = "L6b"
)
dev.off()

# --- Plot 3: Inhibitory subclasses (PVALB, SST, LAMP5, VIP) ---
message("Generating: Inhibitory subclasses...")
svglite(file.path(output_dir, "MTG_inh_subclasses.svg"),
     width = 10, height = 10)
makefigure_human(
  humanTxDb,
  bam1 = pvalb, bam2 = sst, bam3 = lamp5, bam4 = vip,
  chrID = chrID_human, coords = coords_human,
  bam1title = "PVALB", bam2title = "SST",
  bam3title = "LAMP5", bam4title = "VIP"
)
dev.off()

# --- Plot 4: Excitatory vs Inhibitory (simple comparison) ---
message("Generating: Excitatory vs Inhibitory...")
svglite(file.path(output_dir, "MTG_exc_vs_inh.svg"),
     width = 10, height = 6)
makefigure_human(
  humanTxDb,
  bam1 = excitatory, bam2 = inhibitory,
  chrID = chrID_human, coords = coords_human,
  bam1title = "Excitatory", bam2title = "Inhibitory"
)
dev.off()

message("All figures saved to: ", output_dir)
