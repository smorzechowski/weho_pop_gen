# Calculate the cross entropy criterion to pick the best value of K ancestral populations

setwd('/n/netscratch/edwards_lab/Lab/smorzechowski/meliphagid/analysis/2024-11-03/19-lea')


library(LEA)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(stringr)
library(cowplot)
library(gridExtra)



gen.imp <- read.geno("Nleu_autos_lea_depth_filt_imputed_thin_plink.geno")

dim(gen.imp)

################
# running snmf #
################

project.snmf = snmf("Nleu_autos_lea_depth_filt_imputed_thin_plink.geno",
                    K = 1:10, 
                    entropy = TRUE, 
                    repetitions = 10,
                    project = "new")

# using project = "new" saves it to a file with .snmfProject extension

png("sNMF_K1_to_10.png", width = 6, height = 4, units = "in", res = 350)
# plot cross-entropy criterion of all runs of the project
plot(project.snmf, cex = 1.2, col = "lightblue", pch = 19)
dev.off()



