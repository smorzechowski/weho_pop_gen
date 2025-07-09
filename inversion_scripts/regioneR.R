# Enrichment tests
# Sophie M Orzechowski
# March 2025

#install.packages("regioneR")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("regioneR")
#BiocManager::install("rtracklayer")

library(regioneR)
library(GenomicRanges)  # For creating GRanges objects
library(rtracklayer)
library(ggplot2)
library(patchwork)
library(cowplot)

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

data <- data.frame(pt1$numOverlaps$permuted)

pt1

plot1 <- ggplot(data,aes(pt1.numOverlaps.permuted))+
  geom_density(trim=TRUE)+
  xlab("Overlap of avian candidate climate genes in MORs")+
  ylab("Density")+
  ggtitle("A.")+
#  geom_vline(xintercept=79.337,linetype="dashed",color='black')+
  geom_vline(xintercept=98,color='red',linewidth=2)+
  theme_minimal()+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=14),
        title=element_text(size=20))
  

plot1



pt2 <- overlapPermTest(
  A = gr_vertebrate,         # Candidate regions
  B = gr_PCR,          # Annotation regions (e.g., inversions)
  genome = genome_gr, # Specify the genome (or you can provide a custom genome definition)
  ntimes = 1000             # Number of permutations to build the null distribution
)

head(pt2$numOverlaps)


null_mean <- mean(pt2$numOverlaps$permuted)
fold_enrichment <- 143 / null_mean

pt2
fold_enrichment

data <- data.frame(pt2$numOverlaps$permuted)

plot2 <- ggplot(data,aes(pt2.numOverlaps.permuted))+
  geom_density(trim=TRUE)+
  xlab("Overlap of vertebrate candidate climate genes in MORs")+
  ylab("Density")+
  ggtitle("B.")+
 # geom_vline(xintercept=95.497,linetype="dashed",color='black')+
  geom_vline(xintercept=143,color='red',linewidth=2)+
  theme_minimal()+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=14),
        title=element_text(size=20))

final_plot <- plot1/plot2
final_plot


# Dummy plot to create legend
dummy_plot <- ggplot() +
  geom_line(aes(x = Inf, y = Inf, color = "null distribution"), linewidth = 1) +
  geom_line(aes(x = Inf, y = Inf, color = "observed overlap"), linewidth = 1) +
  scale_color_manual(
    name = NULL,
    values = c("null distribution" = "black", "observed overlap" = "red")
  ) +
  theme_void() +
  theme(legend.position = "top",
        legend.text=element_text(size=14))

dummy_plot

# Extract the legend
components <- ggplotGrob(dummy_plot)

# Print all grobs in the plot to see where the legend is
print(components)
custom_legend <- gtable::gtable_filter(components, "guide-box")


final_plot2 <- plot_grid(
  custom_legend,   # Extracted legend
  final_plot,               # Your main ggplot figure
  ncol = 1,        # Arrange them vertically
  rel_heights = c(0.1,1)  # Give more space to the main plot
)

print(final_plot2)


ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/enrichment_avian_vertebrate_candidates.png", 
       plot = final_plot2, dpi = 350, width = 8, height = 8, units = "in")
