# Run RDA
library(vegan)
library(tidyverse)
library(LEA)


gen.imp <- read.geno("Nleu_autos_lea_depth_filt_imputed_thin_plink.geno")
dim(gen.imp)

pred <- read.table("Nleu_bioclim_variables_65ind_autos_ordered_noheader.txt",header=F)
dim(pred)

PC1 <- read.table("Nleu_PC1_65ind_autos_excl_inv.txt",header=F)
colnames(PC1) <- "PC1"

data <- cbind(pred,PC1)

cor(data,method="pearson")


WEHE.rda <- rda(gen.imp ~ V10 + V15 + Condition(PC1), data=data, scale=T)

RsquareAdj(WEHE.rda)

summary(eigenvals(WEHE.rda, model = "constrained"))

screeplot(WEHE.rda)


signif.full <- anova.cca(WEHE.rda, parallel=getOption("mc.cores"))



# Get importance of each environmental variable
#importance <- summary(WEHE.rda)$constr.chi

# Print variable importance
#print(importance)
