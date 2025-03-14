

install.packages("regioneR")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("regioneR")
BiocManager::install("rtracklayer")

library(regioneR)
library(GenomicRanges)  # For creating GRanges objects
library(rtracklayer)

setwd("~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/regioneR")

gr_PCR <- import("inversions_updated.bed",format="BED")

gr_avian <- import("Nleu_avian_candidate_climate_genes_sorted_IDs.bed")
gr_vertebrate <- import("Nleu_vertebrate_pEAGs_genes_sorted_IDs.bed")

genome <- read.table("Nleucotis_hifi_v1.0_hc_sm_fx_scaffolded_PAR_masked.fasta.genome",header=F,comment.char = "")
genome$V2 <- as.numeric(genome$V2)

genome_vector <- setNames(genome$V2,genome$V1)
str(genome_vector)
genome_gr <- GRanges(
  seqnames = names(genome_vector),
  ranges = IRanges(start = 1, end = as.numeric(genome_vector))
)



pt1 <- overlapPermTest(
  A = gr_avian,         # Candidate regions
  B = gr_PCR,          # Annotation regions (e.g., inversions)
  genome = genome_gr, # Specify the genome (or you can provide a custom genome definition)
  ntimes = 1000             # Number of permutations to build the null distribution
)

null_mean <- mean(pt1$numOverlaps$permuted)
fold_enrichment <- 98 / null_mean

fold_enrichment




pt2 <- overlapPermTest(
  A = gr_vertebrate,         # Candidate regions
  B = gr_PCR,          # Annotation regions (e.g., inversions)
  genome = genome_gr, # Specify the genome (or you can provide a custom genome definition)
  ntimes = 1000             # Number of permutations to build the null distribution
)

head(pt$numOverlaps)


null_mean <- mean(pt$numOverlaps$permuted)
fold_enrichment <- 143 / null_mean

fold_enrichment
