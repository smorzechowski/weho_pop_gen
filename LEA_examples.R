# Run LEA
# Examples to understand how LEA works
setwd("~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/LEA")

library(LEA)
library(tidyverse)

# https://rdrr.io/bioc/LEA/man/lfmm2.html

# convert ped file to lfmm
#ped2geno("Nleu_autos_lea_depth_filt_imputed_thin_plink.ped")
#geno2lfmm("Nleu_autos_lea_depth_filt_imputed_thin_plink.geno")

# 221 individuals 128147 loci

# read in lfmm file and environmental predictors
#gen.imp <- read.lfmm("Nleu_autos_lea_depth_filt_imputed_thin_plink.lfmm")
#dim(gen.imp)

############################################################################################
### Example of analysis using lfmm2 ###


# Simulation with 10 target loci, with effect sizes ranging between -10 an 10 
# n = 100 individuals and L = 1000 loci

X <- as.matrix(rnorm(100)) # causal environmental variable
B <- rep(0, 1000) 
target <- sample(1:1000, 10) # target loci
B[target] <- runif(10, -10, +10) # effect sizes

# Creating hidden factors and loadings

U <- t(tcrossprod(as.matrix(c(-1,0.5,1.5)), X))+ matrix(rnorm(300), ncol = 3)
V <- matrix(rnorm(3000), ncol = 3)

# Simulating a binarized matrix containing haploid genotypes 
# Simulation performed with the generative LFMM

Y <- tcrossprod(as.matrix(X), B) + tcrossprod(U, V) + matrix(rnorm(100000, sd = .5), nrow = 100)
Y <- matrix(as.numeric(Y > 0), ncol = 1000)

######################################
# Fitting an LFMM with K = 3 factors #
######################################

mod2 <- lfmm2(input = Y, env = X, K = 3)

# Computing P-values and plotting their minus log10 values 
# Target loci are highlighted

pv <- lfmm2.test(object = mod2, input = Y, env = X, linear = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .4, pch = 19)
points(target, -log10(pv$pvalues[target]), col = "red")


###########################################################################################


# load simulated data
data("offset_example")
# 200 diploid individuals genotyped at 510 SNP 
Y <- offset_example$geno
str(Y)
class(Y)
dim(Y)
# 4 environmental variables
X <- offset_example$env
str(X)
class(X)
dim(X)
mod.lfmm2 <- lfmm2(input = Y, env = X, K = 2)


pv <- lfmm2.test(object = mod.lfmm2,
                 input = Y,
                 env = X,
                 full = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .5, pch = 19)
abline(h = -log10(0.0001), lty = 2, col = "orange")


p.values <- lfmm.pvalues(lfmm.res, causal = TRUE) 
# or use your own approach for multiple testing correction

q.values <- p.adjust(pv$pvalues, method = "BH")
max(q.values)
min(q.values)
max(pv$pvalues)
min(pv$pvalues)
###########################################################################################


library(LEA)
# Creation of a genotype matrix data file: "genotypes.lfmm" # The data include 400 SNPs for 50 individuals.
data("tutorial")

# Write genotypes in the lfmm format
write.lfmm(tutorial.R, "genotypes.lfmm")
# Write genotypes in the geno format

write.geno(tutorial.R, "genotypes.geno")
# creation of an environment gradient file: gradient.env.
# The .env file contains a single ecological variable


################################################################################
#### Example
# load simulated data
data("offset_example")
# 200 diploid individuals genotyped at 510 SNP 
Y <- offset_example$geno
# 4 environmental variables
X <- offset_example$env
str(X)
mod.lfmm2 <- lfmm2(input = Y, env = X, K = 2)


pv <- lfmm2.test(object = mod.lfmm2,
                 input = Y,
                 env = X,
                 full = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .5, pch = 19)
abline(h = -log10(0.1/510), lty = 2, col = "orange")





# Environmental variable X = as.matrix(rnorm(n)) # effect sizes
#B = rep(0, L)