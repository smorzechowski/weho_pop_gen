# Run GEA analysis on White-eared Honeyeater resequencing data
# Sophia MacRae Orzechowski
# January 2025


################################################################################

#gc()
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("LEA")

# convert ped file to geno for lfmm2()
#ped2geno("Nleu_autos_lea_depth_filt_imputed_thin_plink.ped")
#ped2geno("Nleu_autos_neoZ_lea_males_depth_filt_imputed_thin_plink.ped")
#ped2geno("Nleu_autos_lea_depth_filt_imputed_thin_sans_inversions_with_dummy_plink.ped")

# convert text file to env
#pred <- as.matrix(read.table("Nleu_bioclim_variables_65ind_autos_ordered_noheader.txt",header=F))
#pred_matrix <- matrix(pred, nrow = nrow(pred), ncol = ncol(pred))
#write.env(pred_matrix,"Nleu_bioclim_variables_65ind_autos_ordered_noheader.env")

#pred <- as.matrix(read.table("Nleu_bioclim_variables_43males_autos_neoZ_ordered_noheader.txt",header=F))
#pred_matrix <- matrix(pred, nrow = nrow(pred), ncol = ncol(pred))
#write.env(pred_matrix,"Nleu_bioclim_variables_43males_autos_neoZ_ordered_noheader.env")

################################################################################

setwd("/n/netscratch/edwards_lab/Lab/smorzechowski/meliphagid/analysis/2024-11-03/19-lea/")

library(LEA)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(stringr)
library(cowplot)
library(gridExtra)
library(ggrepel)



################################################################################
## All autosomes except new PAR region, 65 individuals


gen.imp <- read.geno("Nleu_autos_lea_depth_filt_imputed_thin_plink.geno")

dim(gen.imp)

#pred <- read.env("Nleu_bioclim_variables_65ind_autos_ordered_noheader.env")
pred <- read.env("Nleu_bioclim_PC1_PC2_65ind_autos.env")


dim(pred)
str(pred)
class(pred) 


map <- read.table("Nleu_autos_lea_depth_filt_imputed_thin_plink.map",header=F,comment.char = "")



################################################################################
# Run LEA models
mod2 <- lfmm2(input = gen.imp, env = pred, K = 2)
#mod3 <- lfmm2(input = gen.imp, env = pred, K = 3)
#mod4 <- lfmm2(input = gen.imp, env = pred, K = 4)




################################################################################
# Computing P-values and plotting their minus log10 values 

# Run multivariate tests, accounting for shared variance (full covariance matrix)
# K=2
pv2.full <- lfmm2.test(object = mod2,
                       input = gen.imp,
                       env = pred,
                       full=TRUE)
# K=3
pv3.full <- lfmm2.test(object = mod3,
                       input = gen.imp,
                       env = pred,
                       full=TRUE)



################################################################################
##### Built-in FDR procedure #####

######
# Calculate adjusted pvalues for the full model
pv2.full.q.values <- p.adjust(pv2.full$pvalues, method = "BH")
pv3.full.q.values <- p.adjust(pv3.full$pvalues, method = "BH")
######



################################################################################
# Target regions across the genome: chromosomes with putative inversions

Chr4_target = seq(from = 474812, to = 583562,by=1)
Chr6_target = seq(from = 583563, to = 650966,by=1)
Chr8_target = seq(from = 700446, to = 757547,by=1)
Chr9_target = seq(from = 757547, to = 773605,by=1)
Chr10_target = seq(from = 139148, to = 158627,by=1)
Chr13_target = seq(from = 172691, to = 201522,by=1)
Chr14_target = seq(from = 201523, to = 220168,by=1)
Chr15_target = seq(from = 220169, to = 230157,by=1)
Chr16_target = seq(from = 230158, to = 247440,by=1)
Chr17_target = seq(from = 247441, to = 263872,by=1)
Chr18_target = seq(from = 263872, to = 272643,by=1)
Con15l_target = seq(from = 773606, to = 789999,by=1)
Con45l_target = seq(from = 790000, to = 807105,by=1)

dim(gen.imp)
dim(pred)


################################################################################
#Create plots
################################################################################
# Full: all variables combined for both K=2 and K=3

# Plot of all variables combined with BH correction K=2
plot(-log10(pv2.full.q.values), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.0000005), lty = 2, col = "darkred")

points(Chr4_target, -log10(pv2.full.q.values[Chr4_target]), col = "red")
points(Chr6_target, -log10(pv2.full.q.values[Chr6_target]), col = "blue")
points(Chr8_target, -log10(pv2.full.q.values[Chr8_target]), col = "lightblue")
points(Chr9_target, -log10(pv2.full.q.values[Chr9_target]), col = "magenta")
points(Chr10_target, -log10(pv2.full.q.values[Chr10_target]), col = "orange")
points(Chr13_target, -log10(pv2.full.q.values[Chr13_target]), col = "green")
points(Chr14_target, -log10(pv2.full.q.values[Chr14_target]), col = "yellow")
points(Chr15_target, -log10(pv2.full.q.values[Chr15_target]), col = "black")
points(Chr16_target, -log10(pv2.full.q.values[Chr16_target]), col = "darkgreen")
points(Chr17_target, -log10(pv2.full.q.values[Chr17_target]), col = "pink")
points(Chr18_target, -log10(pv2.full.q.values[Chr18_target]), col = "brown")
points(Con15l_target, -log10(pv2.full.q.values[Con15l_target]), col = "darkblue")
points(Con45l_target, -log10(pv2.full.q.values[Con45l_target]), col = "darkorange")


# Add a legend to the plot
legend("topleft", 
       legend = c("Chr4", "Chr6", "Chr8", "Chr9", "Chr10", 
                  "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", 
                  "Chr18", "Con15l", "Con45l"),
       col = c("red", "blue", "lightblue", "magenta", "orange", 
               "green", "yellow", "black", "darkgreen", "pink", 
               "brown", "darkblue", "darkorange"),
       pch = 16, # Use circles as point markers
       bty = "n") # Remove box around legend

################################################################################
# Plot of all variables combined with BH correction K=2


# combine map with adjusted p-values
pv2.full.map <- cbind(pv2.full.q.values,map)
table(pv2.full.map$V1[pv2.full.map$pv2.full.q.values<0.0005])
table(pv2.full.map$V1[pv2.full.map$pv2.full.q.values<0.005])

# plot

ggplot(pv2.full.map,
       aes(x = V4/1000000, 
           y = -log10(pv2.full.q.values))) + 
  # 1) Plot ALL points in grey
  geom_point(size = 0.4,# roughly corresponds to cex = 0.4 in base R
             shape = 19)+
  facet_wrap(.~V1,scales="free_x")+
  geom_hline(
    yintercept = -log10(0.0005), 
    linetype = "dashed", 
    color = "darkred"
  )

