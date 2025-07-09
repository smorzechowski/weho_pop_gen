# PCA of LD pruned data from plink

library(ggplot2)
setwd("~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/plink")

#library(tidyverse)
pca <- read.table("Nleu_interval_1_combined_varsite_recode_qd_maf05_filt_relatives_LD.eigenvec",header=F)
eigenval <- read.table("Nleu_interval_1_combined_varsite_recode_qd_maf05_filt_relatives_LD.eigenval",header=F)
# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))


# location
pca$loc <- rep(NA, length(pca$ind))
pca$loc[grep("Ben", pca$ind)] <- "Weddin"
pca$loc[grep("Pill", pca$ind)] <- "Pilliga"
pca$loc[grep("Nom", pca$ind)] <- "Nombinnie"
pca$loc[grep("Reedy", pca$ind)] <- "ReedyCreek"
pca$loc[grep("Walch", pca$ind)] <- "Walcha"
pca$loc[grep("Gund", pca$ind)] <- "Gundagui"
pca$loc[grep("Berht", pca$ind)] <- "Weddin"
pca$loc[grep("Berth", pca$ind)] <- "Weddin"
pca$loc[grep("Inga", pca$ind)] <- "Ingalba"
pca$loc[grep("Binya", pca$ind)] <- "Binya"
pca$loc[grep("Mallee", pca$ind)] <- "Mallee"
pca$loc[grep("Moon", pca$ind)] <- "Moonbi"
pca$loc[grep("Zost", pca$ind)] <- "Zosterops"
pca$loc[grep("Eual", pca$ind)] <- "Weddin"
pca$loc[grep("Bog", pca$ind)] <- "BoggyCreek"
pca$loc[grep("Mull", pca$ind)] <- "Mullions"
pca$loc[grep("Talla", pca$ind)] <- "Tallaganda"

# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

ggplot(pve, aes(PC, V1)) + geom_bar(stat = "identity")+
ylab("Percentage variance explained")+
theme_light()


# calculate the cumulative sum of the percentage variance explained
cumsum(pve$V1)


# plot pca

pca_filt <- pca[pca$loc!="Walcha" & pca$ind!="Nom_100_M"&pca$ind!="Nom_101_F",]
pca_filt <- pca[pca$loc!="Walcha" & pca$loc!="ReedyCreek",]

ggplot(pca_filt, aes(PC1, PC2, col = loc)) + geom_point(size = 3)+
  geom_text(aes(label=ind), vjust = -0.5)
b <- b + scale_colour_manual(values = c("red", "blue"))
b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
