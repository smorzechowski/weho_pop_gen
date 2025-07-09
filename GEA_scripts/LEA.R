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

# Look at specific environmental variables
#BIO3,BIO4, BIO8, BIO15

# "nonredundant"
#pred <- pred[,c(1,2,4,5,6,12,13,14,15)]

# Highest loadings in PC1
#pred <- pred[,c(1,5,12,17)]

# Highest loadings in PC2; cor < 0.7 
#pred <- pred[,c(2,3,8,15)]

# One variable from PC1 and PC2: least correlated
#pred <- pred[,c(10,15)]

# Isothermality
#pred <- pred[,c(3)]

# Mean temperature of warmest quarter
#pred <- pred[,c(10)]


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
# #Evaluate collinearity between latent variables and environment
# Extract the factor matrix U (n x K)
#U <- mod4@U 

# Evaluate collinearity for each environment variable:
#for (j in 1:ncol(pred)) {
#  lm_fit <- lm(pred[,j] ~ U)
#  r2_val <- summary(lm_fit)$r.squared
#  cat(colnames(pred)[j], "-> R^2 with factors:", round(r2_val, 3), "\n")
#}


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

# Run univariate tests for each variable, NOT accounting for shared variance (diagonal covariance matrix) between predictors
# K=2
pv2 <- lfmm2.test(object = mod2,
                       input = gen.imp,
                       env = pred,
                       full=FALSE)
# K=3
pv3 <- lfmm2.test(object = mod3,
                       input = gen.imp,
                       env = pred,
                       full=FALSE)
# K=4
pv4 <- lfmm2.test(object = mod4,
                  input = gen.imp,
                  env = pred,
                  full=FALSE)

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
##### Manual FDR procedure #####
n = 100
#loci
L = 811504
q=0.1 #FDR
level = q * (1) / L 
min(level)
max(level)
w = levelw = which(sort(pvalues) < q * (1:L) / L) 
candidates = order(pvalues)[w]
################################################################################

##### Built-in FDR procedure #####

######
# Calculate adjusted pvalues for the full model
pv2.full.q.values <- p.adjust(pv2.full$pvalues, method = "BH")
pv3.full.q.values <- p.adjust(pv3.full$pvalues, method = "BH")
######




######
# Calculate adjusted p-values for individual variables

# Highest loadings for PC1; K=2
pv2.bio1.q.values <- p.adjust(pv2$pvalues[1,],method="BH")
pv2.bio5.q.values <- p.adjust(pv2$pvalues[2,],method="BH")
pv2.bio12.q.values <- p.adjust(pv2$pvalues[3,],method="BH")
pv2.bio17.q.values <- p.adjust(pv2$pvalues[4,],method="BH")

# Highest loadings for PC1; K=3
pv3.bio1.q.values <- p.adjust(pv3$pvalues[1,],method="BH")
pv3.bio5.q.values <- p.adjust(pv3$pvalues[2,],method="BH")
pv3.bio12.q.values <- p.adjust(pv3$pvalues[3,],method="BH")
pv3.bio17.q.values <- p.adjust(pv3$pvalues[4,],method="BH")

#pred <- pred[,c(2,3,8,15)]
# Highest loadings for PC2; K=2
pv2.bio2.q.values <- p.adjust(pv2$pvalues[1,],method="BH")
pv2.bio3.q.values <- p.adjust(pv2$pvalues[2,],method="BH")
pv2.bio8.q.values <- p.adjust(pv2$pvalues[3,],method="BH")
pv2.bio15.q.values <- p.adjust(pv2$pvalues[4,],method="BH")

# Highest loadings for PC2; K=3
pv3.bio2.q.values <- p.adjust(pv3$pvalues[1,],method="BH")
pv3.bio3.q.values <- p.adjust(pv3$pvalues[2,],method="BH")
pv3.bio8.q.values <- p.adjust(pv3$pvalues[3,],method="BH")
pv3.bio15.q.values <- p.adjust(pv3$pvalues[4,],method="BH")


# One variable each from PC1 and PC2; K=2
pv2.bio10.q.values <- p.adjust(pv2$pvalues[1,],method="BH")
pv2.bio15.q.values <- p.adjust(pv2$pvalues[2,],method="BH")


# One variable each from PC1 and PC2; K=3
pv3.bio10.q.values <- p.adjust(pv3$pvalues[1,],method="BH")
pv3.bio15.q.values <- p.adjust(pv3$pvalues[2,],method="BH")

# One variable each from PC1 and PC2; K=4
pv4.bio10.q.values <- p.adjust(pv4$pvalues[1,],method="BH")
pv4.bio15.q.values <- p.adjust(pv4$pvalues[2,],method="BH")


################################################################################
pv2.bio3.q.values <- p.adjust(pv2$pvalues[1,],method="BH")
pv2.bio2.q.values <- p.adjust(pv2$pvalues[8,],method="BH")
pv2.bio4.q.values <- p.adjust(pv2$pvalues[7,],method="BH")
pv2.bio5.q.values <- p.adjust(pv2$pvalues[6,],method="BH")
pv2.bio6.q.values <- p.adjust(pv2$pvalues[5,],method="BH")
pv2.bio8.q.values <- p.adjust(pv2$pvalues[4,],method="BH")
pv2.bio13.q.values <- p.adjust(pv2$pvalues[3,],method="BH")
pv2.bio14.q.values <- p.adjust(pv2$pvalues[2,],method="BH")




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
## ggplot2 version: Plot of all variables combined with BH correction K=2


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


# Create target_positions variable to color all chromosomes with inversions
pv2.full.map$target_positions <- 0

# Now assign 1 to each specified row-index range
# This colors the entire chromosome containing the inversion
pv2.full.map$target_positions[474812:583562] <- 1
pv2.full.map$target_positions[583563:650966] <- 1
pv2.full.map$target_positions[700446:757547] <- 1
pv2.full.map$target_positions[757547:773605] <- 1
pv2.full.map$target_positions[139148:158627] <- 1
pv2.full.map$target_positions[172691:201522] <- 1
pv2.full.map$target_positions[201523:220168] <- 1
pv2.full.map$target_positions[220169:230157] <- 1
pv2.full.map$target_positions[230158:247440] <- 1
pv2.full.map$target_positions[247441:263872] <- 1
pv2.full.map$target_positions[263872:272643] <- 1
pv2.full.map$target_positions[773606:789999] <- 1
pv2.full.map$target_positions[790000:807105] <- 1


# Create data frame for vertical lines
# Updated breakpoints 26 Feb 2025
lines_df <- data.frame(
  Chr = c(
   # "Chr_4", "Chr_4",
    "Chr_4", "Chr_4",
    "Chr_6", "Chr_6",
   # "Chr_6", "Chr_6",
    "Chr_8", "Chr_8",
    "Chr_8", "Chr_8",
   # "Chr_9", "Chr_9",
    "Chr_9", "Chr_9",
    "pri#ptg000015l", "pri#ptg000015l",
    "Chr_10", "Chr_10",
    "Chr_13", "Chr_13",
    "Chr_14", "Chr_14",
    "Chr_15", "Chr_15",
    "Chr_17", "Chr_17",
    "Chr_18", "Chr_18",
    "pri#ptg000045l", "pri#ptg000045l"
  ),
  vline_pos = c(
    # Chr_4_:95375035-105124963
    # 95375035, 105124963,
    # Chr_4_:107824981-120434650
    # 107824981, 120434650,
    94877058, 120434650, 
    # Chr_6_:8125001-14575038
    #8125001, 14575038,
    # Chr_6_:15924948-328749905
    #15924948, 32874990,
    8125001,32874990,
    # Chr_8_:1-10724985
    1, 10724985,
    # Chr_8_:5597491-37874992
    56196169, 37874992,
    # Chr_9_:16374481-25249745
    #16374481, 2524974,
    # Chr_9_:16374481-18339159
    1631342, 18339159,
    # pri#ptg000015l:6864979-17909547
    6864979, 17909547,
    # Chr_10_:9675047-22345924
    9675047, 22345924,
    # Chr_13_:15724920-25076809
    # Chr_13_:14500000-25076809
    15000000, 23500000,
    # Chr_14_:1-101749905
    1, 10174990,
    # Chr_15_:1-877497
    1, 8774970,
    # Chr_17_:12914902-19617658
    12708285, 19617658,
    # Chr_18_:1865069-10487981
    2865069, 10487981,
    # pri#ptg000045l:1-128348605
    1, 13100247
  ),
  stringsAsFactors = TRUE
)




intervals_df <- lines_df %>%
  group_by(Chr) %>%
  arrange(vline_pos) %>%
  # Create pairs: rows 1-2 form the first interval, rows 3-4 form the second, etc.
  mutate(row_index = row_number()) %>%
  # For each pair, the "start" is an odd row_index, the "end" is the following even row_index
  # We'll pivot those out into one row per pair.
  summarize(
    start = vline_pos[row_index %% 2 == 1], 
    end   = vline_pos[row_index %% 2 == 0],
    .groups = "keep"
  ) %>%
  # It's often handy to reorder columns or rename them
  select(Chr, start, end)


pv2.full.map$Chr <- as.character(pv2.full.map$V1)
intervals_df$Chr <- as.character(intervals_df$Chr)

pv2.full.map$in_interval <- FALSE

pv2.full.map$Chr <- str_replace(pv2.full.map$Chr,'scaffold','Chr')
pv2.full.map$Chr <- str_remove(pv2.full.map$Chr,'_RagTag')

# For each interval, mark those points as in_interval
for(i in seq_len(nrow(intervals_df))) {
  chr_i   <- intervals_df$Chr[i]
  start_i <- intervals_df$start[i]
  end_i   <- intervals_df$end[i]
  
  pv2.full.map$in_interval <- pv2.full.map$in_interval | 
    (pv2.full.map$Chr == chr_i & pv2.full.map$V4 >= start_i & pv2.full.map$V4 <= end_i)
}


table(pv2.full.map$Chr)

pv2.full.map$Chr <- factor(pv2.full.map$Chr,levels=c("Chr_1","Chr_3","Chr_4","Chr_6",
                                    "Chr_7","Chr_8","Chr_9","Chr_10",
                                    "Chr_11","Chr_12","Chr_13","Chr_14",
                                    "Chr_15","Chr_16","Chr_17","Chr_18",
                                    "Chr_19","Chr_20","Chr_21","Chr_22",
                                    "Chr_23","Chr_24","Chr_25","Chr_26",
                                    "Chr_27","Chr_28","Chr_29","Chr_30",
                                    "Chr_31","Chr_32","Chr_33","Chr_34",
                                    "Chr_35","Chr_36","Chr_37","Chr_38",
                                    "Chr_39","pri#ptg000015l","pri#ptg000031l",
                                    "pri#ptg000045l","pri#ptg000050l","pri#ptg000051l",
                                    "pri#ptg000085l","pri#ptg000090l"))


# fix target positions
pv2.full.map$target_positions[pv2.full.map$V1=="pri#ptg000031l"] <- 0
pv2.full.map$target_positions[pv2.full.map$V1=="scaffold_16_RagTag"] <- 0

# position variable for geom_segment()
y_top <- max(-log10(pv2.full.map$pv2.full.q.values), na.rm = TRUE) * 1.1
y_mid <- max(-log10(pv2.full.map$pv2.full.q.values), na.rm = TRUE) * 0.5




# add candidate genes

climategenes <- read.table('candidate_climate_inversion_genes.txt',header=F,strip.white = TRUE)
colnames(climategenes)[1] <- "gene"
#geneloc <- read.table('pv2.full.GEA.hits_inversions_25kb_windows.txt',header=F,comment.char = "",sep='\t')
#geneloc <- read.table('pv2.full.snps_25kb_window_genes.bed',header=F,comment.char = "",sep='\t')
geneloc <- read.delim('Nleu_all_genes.bed',header=F,comment.char = "",sep="\t",strip.white = TRUE,fill=TRUE,quote="")

#geneloc$gene <- str_extract(geneloc$V7, "(?<=ID=)[^_]+")
geneloc$gene <- str_extract(geneloc$V4, "(?<=ID=)[^_]+")
geneloc_sort <- geneloc[order(geneloc$gene),]
geneloc_sort_nondup <- geneloc_sort[!duplicated(geneloc_sort$gene),]
climate_geneloc <- geneloc_sort[is.element(geneloc_sort$gene,climategenes$V1),]


climate_geneloc <- left_join(climategenes,geneloc)
climate_geneloc$Chr <- str_replace(climate_geneloc$V1,'scaffold','Chr')
climate_geneloc$Chr <- str_remove(climate_geneloc$Chr,'_RagTag')

climate_geneloc$Chr <- factor(climate_geneloc$Chr,levels=c("Chr_1","Chr_3","Chr_4","Chr_6",
                                                     "Chr_7","Chr_8","Chr_9","Chr_10",
                                                     "Chr_11","Chr_12","Chr_13","Chr_14",
                                                     "Chr_15","Chr_16","Chr_17","Chr_18",
                                                     "Chr_19","Chr_20","Chr_21","Chr_22",
                                                     "Chr_23","Chr_24","Chr_25","Chr_26",
                                                     "Chr_27","Chr_28","Chr_29","Chr_30",
                                                     "Chr_31","Chr_32","Chr_33","Chr_34",
                                                     "Chr_35","Chr_36","Chr_37","Chr_38",
                                                     "Chr_39","pri#ptg000015l","pri#ptg000031l",
                                                     "pri#ptg000045l","pri#ptg000050l","pri#ptg000051l",
                                                     "pri#ptg000085l","pri#ptg000090l"))



# plot ggplot with inversions displayed as blue bars

p <- ggplot(pv2.full.map,
       aes(x = V4/1000000, 
           y = -log10(pv2.full.q.values),color=in_interval)) + 
  # 1) Plot ALL points in grey
  geom_point(size = 0.4,# roughly corresponds to cex = 0.4 in base R
    shape = 19
  ) +
  geom_segment(data=intervals_df,aes(x=start/1e6,xend=end/1e6),y=y_top,yend=y_top,color='blue',linewidth=2,inherit.aes=FALSE)+
  geom_point(data=climate_geneloc,aes(x=V2/1e6,y=y_mid),color='red',inherit.aes=FALSE)+
#  geom_vline(
#    data = lines_df,
#    aes(xintercept = vline_pos / 1e6),
#    color = "blue",
#    linetype = "dashed",
#  )+
  facet_wrap(.~Chr,scales="free_x")+
  # 2) Add the horizontal threshold line
  geom_hline(
    yintercept = -log10(0.005), 
    linetype = "dashed", 
    color = "darkred"
  ) +
  # 4) Add a nice theme, axes labels, etc.
  theme_minimal() +
  labs(
    x = "Chromosome position (Mb)",
    y = expression(-log[10](italic(p)))
  )+
  scale_color_manual(values=c("grey","darkblue"))+
  theme(legend.position = "none",
        strip.text = element_text(size=13))

p

# Save the plot as a PNG with 400 dpi, for example 10 x 7 inches:
ggsave("pv2.full.map_PC1_2.png", plot = p, dpi = 400, width = 10, height = 7.5, units = "in")









################################################################################
# Create a list of all SNPs and their locations with adjusted p-value <0.005 (to be moderately conservative)


# combine map with adjusted p-values
pv2.full.map <- cbind(pv2.full.q.values,map)

sigSNPs <- pv2.full.map[pv2.full.map$pv2.full.q.values<0.005,]


write.table(sigSNPs,file="pv2.full.map.sigSNPs_0.005_threshold.txt",quote=FALSE,row.names=FALSE,sep="\t")



################################################################################
# Plot of all variables combined with BH correction K=3
plot(-log10(pv3.full.q.values), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.0000005), lty = 2, col = "darkred")


points(Chr4_target, -log10(pv3.full.q.values[Chr4_target]), col = "red")
points(Chr6_target, -log10(pv3.full.q.values[Chr6_target]), col = "blue")
points(Chr8_target, -log10(pv3.full.q.values[Chr8_target]), col = "lightblue")
points(Chr9_target, -log10(pv3.full.q.values[Chr9_target]), col = "magenta")
points(Chr10_target, -log10(pv3.full.q.values[Chr10_target]), col = "orange")
points(Chr13_target, -log10(pv3.full.q.values[Chr13_target]), col = "green")
points(Chr14_target, -log10(pv3.full.q.values[Chr14_target]), col = "yellow")
points(Chr15_target, -log10(pv3.full.q.values[Chr15_target]), col = "black")
points(Chr16_target, -log10(pv3.full.q.values[Chr16_target]), col = "darkgreen")
points(Chr17_target, -log10(pv3.full.q.values[Chr17_target]), col = "pink")
points(Chr18_target, -log10(pv3.full.q.values[Chr18_target]), col = "brown")
points(Con15l_target, -log10(pv3.full.q.values[Con15l_target]), col = "darkblue")
points(Con45l_target, -log10(pv3.full.q.values[Con45l_target]), col = "darkorange")


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
# PC axes for both K=2 and K=3

#K=2
pv2.PC1.q.values <- p.adjust(pv2$pvalues[1,],method="BH")
plot(-log10(pv2.PC1.q.values), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.00005), lty = 2, col = "darkred")

pv2.PC2.q.values <- p.adjust(pv2$pvalues[2,],method="BH")
plot(-log10(pv2.PC2.q.values), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.00005), lty = 2, col = "darkred")


#K=3
pv3.PC1.q.values <- p.adjust(pv3$pvalues[1,],method="BH")
plot(-log10(pv3.PC1.q.values), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.00005), lty = 2, col = "darkred")

pv3.PC2.q.values <- p.adjust(pv3$pvalues[2,],method="BH")
plot(-log10(pv3.PC2.q.values), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.00005), lty = 2, col = "darkred")



################################################################################
# Highest loading PC1 K=2

# Plot of one variable: bio 1
plot(-log10(pv2.bio1.q.values), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.00005), lty = 2, col = "darkred")

# Plot of one variable: bio 5
plot(-log10(pv2.bio5.q.values), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.00005), lty = 2, col = "darkred")

# Plot of one variable: bio 12
plot(-log10(pv2.bio12.q.values), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.00005), lty = 2, col = "darkred")

# Plot of one variable: bio 15
plot(-log10(pv2.bio17.q.values), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.00005), lty = 2, col = "darkred")


################################################################################
# Highest loading PC2 K=2

# Plot of one variable: bio 2
plot(-log10(pv2.bio2.q.values), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.00005), lty = 2, col = "darkred")

# Plot of one variable: bio 3
plot(-log10(pv2.bio3.q.values), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.00005), lty = 2, col = "darkred")

# Plot of one variable: bio 8
plot(-log10(pv2.bio8.q.values), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.00005), lty = 2, col = "darkred")

# Plot of one variable: bio 15
plot(-log10(pv2.bio15.q.values), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.00005), lty = 2, col = "darkred")

################################################################################
# One variable each from PC1 and PC2 K=2

# Plot of one variable from PC1
plot(-log10(pv2.bio10.q.values), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.00005), lty = 2, col = "darkred")

# Plot of one variable from PC2
plot(-log10(pv2.bio15.q.values), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.00005), lty = 2, col = "darkred")


# One variable each from PC1 and PC2 K=3

# Plot of one variable from PC1
plot(-log10(pv3.bio10.q.values), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.00005), lty = 2, col = "darkred")

points(Chr4_target, -log10(pv3.bio10.q.values[Chr4_target]), col = "red")
points(Chr6_target, -log10(pv3.bio10.q.values[Chr6_target]), col = "blue")
points(Chr8_target, -log10(pv3.bio10.q.values[Chr8_target]), col = "lightblue")
points(Chr9_target, -log10(pv3.bio10.q.values[Chr9_target]), col = "magenta")
points(Chr10_target, -log10(pv3.bio10.q.values[Chr10_target]), col = "orange")
points(Chr13_target, -log10(pv3.bio10.q.values[Chr13_target]), col = "green")
points(Chr14_target, -log10(pv3.bio10.q.values[Chr14_target]), col = "yellow")
points(Chr15_target, -log10(pv3.bio10.q.values[Chr15_target]), col = "black")
points(Chr16_target, -log10(pv3.bio10.q.values[Chr16_target]), col = "darkgreen")
points(Chr17_target, -log10(pv3.bio10.q.values[Chr17_target]), col = "pink")
points(Chr18_target, -log10(pv3.bio10.q.values[Chr18_target]), col = "brown")
points(Con15l_target, -log10(pv3.bio10.q.values[Con15l_target]), col = "darkblue")
points(Con45l_target, -log10(pv3.bio10.q.values[Con45l_target]), col = "darkorange")




# Plot of one variable from PC2
plot(-log10(pv3.bio15.q.values), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.00005), lty = 2, col = "darkred")


# One variable each from PC1 and PC2 K=4

# Plot of one variable from PC1
plot(-log10(pv4.bio10.q.values), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.00005), lty = 2, col = "darkred")

# Plot of one variable from PC2
plot(-log10(pv4.bio15.q.values), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.00005), lty = 2, col = "darkred")

points(Chr4_target, -log10(pv4.bio15.q.values[Chr4_target]), col = "red")
points(Chr6_target, -log10(pv4.bio15.q.values[Chr6_target]), col = "blue")
points(Chr8_target, -log10(pv4.bio15.q.values[Chr8_target]), col = "lightblue")
points(Chr9_target, -log10(pv4.bio15.q.values[Chr9_target]), col = "magenta")
points(Chr10_target, -log10(pv4.bio15.q.values[Chr10_target]), col = "orange")
points(Chr13_target, -log10(pv4.bio15.q.values[Chr13_target]), col = "green")
points(Chr14_target, -log10(pv4.bio15.q.values[Chr14_target]), col = "yellow")
points(Chr15_target, -log10(pv4.bio15.q.values[Chr15_target]), col = "black")
points(Chr16_target, -log10(pv4.bio15.q.values[Chr16_target]), col = "darkgreen")
points(Chr17_target, -log10(pv4.bio15.q.values[Chr17_target]), col = "pink")
points(Chr18_target, -log10(pv4.bio15.q.values[Chr18_target]), col = "brown")
points(Con15l_target, -log10(pv4.bio15.q.values[Con15l_target]), col = "darkblue")
points(Con45l_target, -log10(pv4.bio15.q.values[Con45l_target]), col = "darkorange")




################################################################################


# Plot of one variable: bio 4
plot(-log10(pv$pvalues[1,]), col = "grey", cex = .4, pch = 19)
pv2.bio4.q.values <- p.adjust(pv2$pvalues[1,],method="BH")
plot(-log10(pv2.bio4.q.values), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.001), lty = 2, col = "darkred")

# Plot of one variable: bio 12
plot(-log10(pv$pvalues[2,]), col = "grey", cex = .4, pch = 19)
pv2.bio12.q.values <- p.adjust(pv2$pvalues[2,],method="BH")
plot(-log10(pv2.bio12.q.values), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.001), lty = 2, col = "darkred")


#K=3
# Plot of one variable: bio 4
plot(-log10(pv3$pvalues[1,]), col = "grey", cex = .4, pch = 19)
pv3.bio4.q.values <- p.adjust(pv3$pvalues[1,],method="BH")
plot(-log10(pv3.bio4.q.values), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.001), lty = 2, col = "darkred")

# Plot of one variable: bio 12
plot(-log10(pv3$pvalues[2,]), col = "grey", cex = .4, pch = 19)
pv3.bio12.q.values <- p.adjust(pv3$pvalues[2,],method="BH")
plot(-log10(pv3.bio12.q.values), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.001), lty = 2, col = "darkred")



################################################################################
# Look at environmental variables separately with adjusted p-values

# Plot of one variable: bio 15
plot(-log10(pv2$pvalues[1,]), col = "grey", cex = .4, pch = 19)
plot(-log10(pv2.bio15.q.values), col = "grey", cex = .4, pch = 19)


# Plot of one variable: bio 14
plot(-log10(pv2$pvalues[2,]), col = "grey", cex = .4, pch = 19)
plot(-log10(pv2.bio14.q.values), col = "grey", cex = .4, pch = 19)

# Plot of one variable: bio 13
plot(-log10(pv2$pvalues[3,]), col = "grey", cex = .4, pch = 19)
plot(-log10(pv2.bio13.q.values), col = "grey", cex = .4, pch = 19)


# Plot of one variable: bio 6
plot(-log10(pv2$pvalues[5,]), col = "grey", cex = .4, pch = 19)
plot(-log10(pv2.bio6.q.values), col = "grey", cex = .4, pch = 19)


# Plot of one variable: bio 5
plot(-log10(pv2$pvalues[6,]), col = "grey", cex = .4, pch = 19)
plot(-log10(pv2.bio5.q.values), col = "grey", cex = .4, pch = 19)


# Plot of one variable: bio 4
plot(-log10(pv2$pvalues[7,]), col = "grey", cex = .4, pch = 19)
plot(-log10(pv2.bio4.q.values), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.0000005), lty = 2, col = "darkred")


# Plot of one variable: bio 2
plot(-log10(pv2$pvalues[8,]), col = "grey", cex = .4, pch = 19)
plot(-log10(pv2.bio2.q.values), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.00005), lty = 2, col = "darkred")


# Plot of one variable: bio 1
plot(-log10(pv2$pvalues[9,]), col = "grey", cex = .4, pch = 19)
plot(-log10(pv2.bio1.q.values), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.00005), lty = 2, col = "darkred")

points(Chr4_target, -log10(pv2.bio1.q.values[Chr4_target]), col = "red")
points(Chr6_target, -log10(pv2.bio1.q.values[Chr6_target]), col = "blue")
points(Chr8_target, -log10(pv2.bio1.q.values[Chr8_target]), col = "lightblue")
points(Chr9_target, -log10(pv2.bio1.q.values[Chr9_target]), col = "magenta")

points(Chr10_target, -log10(pv2.bio1.q.values[Chr10_target]), col = "orange")
points(Chr13_target, -log10(pv2.bio1.q.values[Chr13_target]), col = "green")
points(Chr14_target, -log10(pv2.bio1.q.values[Chr14_target]), col = "yellow")
points(Chr15_target, -log10(pv2.bio1.q.values[Chr15_target]), col = "black")
points(Chr16_target, -log10(pv2.bio1.q.values[Chr16_target]), col = "darkgreen")
points(Chr17_target, -log10(pv2.bio1.q.values[Chr17_target]), col = "pink")
points(Chr18_target, -log10(pv2.bio1.q.values[Chr18_target]), col = "brown")

points(Con15l_target, -log10(pv2.bio1.q.values[Con15l_target]), col = "darkblue")
points(Con45l_target, -log10(pv2.bio1.q.values[Con45l_target]), col = "darkorange")


plot(-log10(pv.mod3$pvalues), col = "grey", cex = .4, pch = 19)

# Bonferonni corrections -check into this more 
abline(h = -log10(0.1/811504), lty = 2, col = "orange")
abline(h = -log10(0.0005/811504), lty = 2, col = "red")
abline(h = -log10(0.0000005/811504), lty = 2, col = "darkred")


points(Chr4_target, -log10(pv$pvalues[Chr4_target]), col = "red")
points(Chr6_target, -log10(pv$pvalues[Chr6_target]), col = "blue")
points(Chr8_target, -log10(pv$pvalues[Chr8_target]), col = "lightblue")
points(Chr9_target, -log10(pv$pvalues[Chr9_target]), col = "magenta")

points(Chr10_target, -log10(pv$pvalues[Chr10_target]), col = "orange")
points(Chr13_target, -log10(pv$pvalues[Chr13_target]), col = "green")
points(Chr14_target, -log10(pv$pvalues[Chr14_target]), col = "yellow")
points(Chr15_target, -log10(pv$pvalues[Chr15_target]), col = "black")
points(Chr16_target, -log10(pv$pvalues[Chr16_target]), col = "darkgreen")
points(Chr17_target, -log10(pv$pvalues[Chr17_target]), col = "pink")
points(Chr18_target, -log10(pv$pvalues[Chr18_target]), col = "brown")

points(Con15l_target, -log10(pv$pvalues[Con15l_target]), col = "darkblue")
points(Con45l_target, -log10(pv$pvalues[Con45l_target]), col = "darkorange")



################################################################################
## All autosomes and neo-Z, including new PAR, 43 males
################################################################################



Mgen.imp <- read.geno("Nleu_autos_neoZ_lea_males_depth_filt_imputed_thin_plink.geno")
dim(Mgen.imp)


#Mpred <- read.env("Nleu_bioclim_variables_43males_autos_neoZ_ordered_noheader.env")
Mpred <- read.env("Nleu_bioclim_PC1_PC2_43males_autos_neoZ.env")

# Isothermality
#Mpred <- Mpred[,c(3)]


# Mean temperature of warmest quarter
#Mpred <- Mpred[,c(10)]


#Mpred <- Mpred[,c(19,8)]

dim(Mgen.imp)
dim(Mpred)
str(Mpred)
class(Mpred) 

Mmap <- read.table("Nleu_autos_neoZ_lea_males_depth_filt_imputed_thin_plink.map",header=F,comment.char = "")



################################################################################
# Run LEA models

Mmod2 <- lfmm2(input = Mgen.imp, env = Mpred, K = 2)
#Mmod3 <- lfmm2(input = Mgen.imp, env = Mpred, K = 3)

################################################################################
# evaluate collinearity between latent factors and environmental variables
# Extract the factor matrix U (n x K)

#U <- Mmod2@U 

# Evaluate collinearity for each environment variable:
#for (j in 1:ncol(Mpred)) {
#  lm_fit <- lm(Mpred[,j] ~ U)
#  r2_val <- summary(lm_fit)$r.squared
#  cat(colnames(Mpred)[j], "-> R^2 with factors:", round(r2_val, 3), "\n")
#}


################################################################################
# Computing P-values and plotting their minus log10 values 
# Target loci are highlighted


# Multivariate tests
Mpv2.full <- lfmm2.test(object = Mmod2,
                        input = Mgen.imp,
                        env = Mpred,
                        full=TRUE)

Mpv3.full <- lfmm2.test(object = Mmod3,
                        input = Mgen.imp,
                        env = Mpred,
                        full=TRUE)

# Univariate tests
Mpv2 <- lfmm2.test(object = Mmod2,
                 input = Mgen.imp,
                 env = Mpred,
                 full=FALSE)

Mpv3 <- lfmm2.test(object = Mmod3,
                      input = Mgen.imp,
                      env = Mpred,
                      full=FALSE)


#individuals
n = 43
#loci
L = 927011

################################################################################
# Adjust p-values using BH method

Mpv2.full.q.values <- p.adjust(Mpv2.full$pvalues,method="BH")
Mpv3.full.q.values <- p.adjust(Mpv3.full$pvalues,method="BH")



################################################################################
# Target regions across the genome: chromosomes with putative inversions
ChrZ_target = seq(from = 139285, to = 253200,by=1)
Chr4_target = seq(from = 589810, to = 698719,by=1)
Chr6_target = seq(from = 698720, to = 766209,by=1)
Chr8_target = seq(from = 815723, to = 872891,by=1)
Chr9_target = seq(from = 872892, to = 888966,by=1)
Chr10_target = seq(from = 253200, to = 272686,by=1)
Chr13_target = seq(from = 286953, to = 315855,by=1)
Chr14_target = seq(from = 315856, to = 334550,by=1)
Chr15_target = seq(from = 334551, to = 344541,by=1)
Chr16_target = seq(from = 344542, to = 361872,by=1)
Chr17_target = seq(from = 361873, to = 378350,by=1)
Chr18_target = seq(from = 378351, to = 387148,by=1)
Con15l_target = seq(from = 888967, to = 905053,by=1)
Con45l_target = seq(from = 905434, to = 922591,by=1)




######################################################################
#Create plots
######################################################################

plot(-log10(Mpv$pvalues), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.00005/927011), lty = 2, col = "darkred")
points(ChrZ_target, -log10(Mpv$pvalues[ChrZ_target]), col = "darkred")
points(Chr4_target, -log10(Mpv$pvalues[Chr4_target]), col = "red")
points(Chr6_target, -log10(Mpv$pvalues[Chr6_target]), col = "blue")
points(Chr8_target, -log10(Mpv$pvalues[Chr8_target]), col = "lightblue")
points(Chr9_target, -log10(Mpv$pvalues[Chr9_target]), col = "magenta")
points(Chr10_target, -log10(Mpv$pvalues[Chr10_target]), col = "orange")
points(Chr13_target, -log10(Mpv$pvalues[Chr13_target]), col = "green")
points(Chr14_target, -log10(Mpv$pvalues[Chr14_target]), col = "yellow")
points(Chr15_target, -log10(Mpv$pvalues[Chr15_target]), col = "black")
points(Chr16_target, -log10(Mpv$pvalues[Chr16_target]), col = "darkgreen")
points(Chr17_target, -log10(Mpv$pvalues[Chr17_target]), col = "pink")
points(Chr18_target, -log10(Mpv$pvalues[Chr18_target]), col = "brown")
points(Con15l_target, -log10(Mpv$pvalues[Con15l_target]), col = "darkblue")
points(Con45l_target, -log10(Mpv$pvalues[Con45l_target]), col = "darkorange")


# Add a legend with filled points
legend("topright", 
       legend = c("ChrZ", "Chr4", "Chr6", "Chr8", "Chr9", "Chr10", 
                  "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", 
                  "Chr18", "Con15l", "Con45l"),
       col = c("darkred", "red", "blue", "lightblue", "magenta", "orange", 
               "green", "yellow", "black", "darkgreen", "pink", 
               "brown", "darkblue", "darkorange"),
       pch = 16,  # Use filled circles
       bty = "n") # Remove box around legend



########################################################################
# PC axes of variation

# full = TRUE, K=2
Mpv2.full.q.values <- p.adjust(Mpv2.full$pvalues,method="BH")
plot(-log10(Mpv2.full.q.values), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.00005), lty = 2, col = "darkred")
points(ChrZ_target, -log10(Mpv2.full.q.values[ChrZ_target]), col = "darkred")
points(Chr4_target, -log10(Mpv2.full.q.values[Chr4_target]), col = "red")
points(Chr6_target, -log10(Mpv2.full.q.values[Chr6_target]), col = "blue")
points(Chr8_target, -log10(Mpv2.full.q.values[Chr8_target]), col = "lightblue")
points(Chr9_target, -log10(Mpv2.full.q.values[Chr9_target]), col = "magenta")
points(Chr10_target, -log10(Mpv2.full.q.values[Chr10_target]), col = "orange")
points(Chr13_target, -log10(Mpv2.full.q.values[Chr13_target]), col = "green")
points(Chr14_target, -log10(Mpv2.full.q.values[Chr14_target]), col = "yellow")
points(Chr15_target, -log10(Mpv2.full.q.values[Chr15_target]), col = "black")
points(Chr16_target, -log10(Mpv2.full.q.values[Chr16_target]), col = "darkgreen")
points(Chr17_target, -log10(Mpv2.full.q.values[Chr17_target]), col = "pink")
points(Chr18_target, -log10(Mpv2.full.q.values[Chr18_target]), col = "brown")
points(Con15l_target, -log10(Mpv2.full.q.values[Con15l_target]), col = "darkblue")
points(Con45l_target, -log10(Mpv2.full.q.values[Con45l_target]), col = "darkorange")


# Add a legend with filled points
legend("topleft", 
       legend = c("ChrZ", "Chr4", "Chr6", "Chr8", "Chr9", "Chr10", 
                  "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", 
                  "Chr18", "Con15l", "Con45l"),
       col = c("darkred", "red", "blue", "lightblue", "magenta", "orange", 
               "green", "yellow", "black", "darkgreen", "pink", 
               "brown", "darkblue", "darkorange"),
       pch = 16,  # Use filled circles
       bty = "n") # Remove box around legend


################################################################################
# ggplot2 version: Plot of all variables combined with BH correction K=2


Mpv2.full.map <- cbind(Mpv2.full.q.values,Mmap)
table(Mpv2.full.map$V1[Mpv2.full.map$Mpv2.full.q.values<0.005])

neoZ_results <- Mpv2.full.map[Mpv2.full.map$V1=="scaffold_2_RagTag" & Mpv2.full.map$Mpv2.full.q.values<0.005,]

# Updated SV breakpoints 26 Feb 2025
lines_df <- data.frame(
  Chr = c(
    # "Chr_4", "Chr_4",
    "Chr_4", "Chr_4",
    "Chr_6", "Chr_6",
    # "Chr_6", "Chr_6",
    "Chr_8", "Chr_8",
    "Chr_8", "Chr_8",
    # "Chr_9", "Chr_9",
    "Chr_9", "Chr_9",
    "pri#ptg000015l", "pri#ptg000015l",
    "Chr_10", "Chr_10",
    "Chr_13", "Chr_13",
    "Chr_14", "Chr_14",
    "Chr_15", "Chr_15",
    "Chr_17", "Chr_17",
    "Chr_18", "Chr_18",
    "pri#ptg000045l", "pri#ptg000045l"
  ),
  vline_pos = c(
    # Chr_4_:95375035-105124963
    # 95375035, 105124963,
    # Chr_4_:107824981-120434650
    # 107824981, 120434650,
    94877058, 120434650, 
    # Chr_6_:8125001-14575038
    #8125001, 14575038,
    # Chr_6_:15924948-328749905
    #15924948, 32874990,
    8125001,32874990,
    # Chr_8_:1-10724985
    1, 10724985,
    # Chr_8_:5597491-37874992
    56196169, 37874992,
    # Chr_9_:16374481-25249745
    #16374481, 2524974,
    # Chr_9_:16374481-18339159
    1631342, 18339159,
    # pri#ptg000015l:6864979-17909547
    6864979, 17909547,
    # Chr_10_:9675047-22345924
    9675047, 22345924,
    # Chr_13_:15724920-25076809
    # Chr_13_:14500000-25076809
    15000000, 23500000,
    # Chr_14_:1-101749905
    1, 10174990,
    # Chr_15_:1-877497
    1, 8774970,
    # Chr_17_:12914902-19617658
    12708285, 19617658,
    # Chr_18_:1865069-10487981
    2865069, 10487981,
    # pri#ptg000045l:1-128348605
    1, 13100247
  ),
  stringsAsFactors = TRUE
)


Mpv2.full.map$Chr <- as.character(Mpv2.full.map$V1)
lines_df$Chr <- as.character(lines_df$Chr)

intervals_df <- lines_df %>%
  group_by(Chr) %>%
  arrange(vline_pos) %>%
  # Create pairs: rows 1-2 form the first interval, rows 3-4 form the second, etc.
  mutate(row_index = row_number()) %>%
  # For each pair, the "start" is an odd row_index, the "end" is the following even row_index
  # We'll pivot those out into one row per pair.
  summarize(
    start = vline_pos[row_index %% 2 == 1], 
    end   = vline_pos[row_index %% 2 == 0],
    .groups = "keep"
  ) %>%
  # It's often handy to reorder columns or rename them
  select(Chr, start, end)



Mpv2.full.map$in_interval <- FALSE

Mpv2.full.map$Chr <- str_replace(Mpv2.full.map$Chr,'scaffold','Chr')
Mpv2.full.map$Chr <- str_remove(Mpv2.full.map$Chr,'_RagTag')

# For each interval, mark those points as in_interval
for(i in seq_len(nrow(intervals_df))) {
  chr_i   <- intervals_df$Chr[i]
  start_i <- intervals_df$start[i]
  end_i   <- intervals_df$end[i]
  
  Mpv2.full.map$in_interval <- Mpv2.full.map$in_interval | 
    (Mpv2.full.map$Chr == chr_i & Mpv2.full.map$V4 >= start_i & Mpv2.full.map$V4 <= end_i)
}

Mpv2.full.map$Chr[Mpv2.full.map$V1=="scaffold_2_RagTag"]<-"Neo_Z"

table(Mpv2.full.map$Chr)

Mpv2.full.map$Chr <- factor(Mpv2.full.map$Chr,levels=c("Neo_Z","Chr_1","Chr_3","Chr_4","Chr_6",
                                                     "Chr_7","Chr_8","Chr_9","Chr_10",
                                                     "Chr_11","Chr_12","Chr_13","Chr_14",
                                                     "Chr_15","Chr_16","Chr_17","Chr_18",
                                                     "Chr_19","Chr_20","Chr_21","Chr_22",
                                                     "Chr_23","Chr_24","Chr_25","Chr_26",
                                                     "Chr_27","Chr_28","Chr_29","Chr_30",
                                                     "Chr_31","Chr_32","Chr_33","Chr_34",
                                                     "Chr_35","Chr_36","Chr_37","Chr_38",
                                                     "Chr_39","pri#ptg000015l","pri#ptg000031l",
                                                     "pri#ptg000045l","pri#ptg000050l","pri#ptg000051l",
                                                     "pri#ptg000085l","pri#ptg000090l"))






lines_df$Chr <- factor(lines_df$Chr,levels=c("Neo_Z","Chr_1","Chr_3","Chr_4","Chr_6",
                                             "Chr_7","Chr_8","Chr_9","Chr_10",
                                             "Chr_11","Chr_12","Chr_13","Chr_14",
                                             "Chr_15","Chr_16","Chr_17","Chr_18",
                                             "Chr_19","Chr_20","Chr_21","Chr_22",
                                             "Chr_23","Chr_24","Chr_25","Chr_26",
                                             "Chr_27","Chr_28","Chr_29","Chr_30",
                                             "Chr_31","Chr_32","Chr_33","Chr_34",
                                             "Chr_35","Chr_36","Chr_37","Chr_38",
                                             "Chr_39","pri#ptg000015l","pri#ptg000031l",
                                             "pri#ptg000045l","pri#ptg000050l","pri#ptg000051l",
                                             "pri#ptg000085l","pri#ptg000090l"))


intervals_df$Chr <- factor(intervals_df$Chr,levels=c("Neo_Z","Chr_1","Chr_3","Chr_4","Chr_6",
                                                     "Chr_7","Chr_8","Chr_9","Chr_10",
                                                     "Chr_11","Chr_12","Chr_13","Chr_14",
                                                     "Chr_15","Chr_16","Chr_17","Chr_18",
                                                     "Chr_19","Chr_20","Chr_21","Chr_22",
                                                     "Chr_23","Chr_24","Chr_25","Chr_26",
                                                     "Chr_27","Chr_28","Chr_29","Chr_30",
                                                     "Chr_31","Chr_32","Chr_33","Chr_34",
                                                     "Chr_35","Chr_36","Chr_37","Chr_38",
                                                     "Chr_39","pri#ptg000015l","pri#ptg000031l",
                                                     "pri#ptg000045l","pri#ptg000050l","pri#ptg000051l",
                                                     "pri#ptg000085l","pri#ptg000090l"))








# plot

# Create a categorical variable for SNPs that are significant
Mpv2.full.map$pcat <- FALSE
Mpv2.full.map$pcat[Mpv2.full.map$Mpv2.full.q.values<0.005] <- TRUE


# position variable for geom_segment()
y_top <- max(-log10(Mpv2.full.map$Mpv2.full.q.values), na.rm = TRUE) * 1.01


# Plot intervals for neo-Z
intervals_neoZ <- data.frame(Chr=c("Neo_Z","Neo_Z","Neo_Z"),start=c(0,74.4,111.6),end=c(74.4,111.6,128.9))
intervals_neoZ$Chr <- as.factor(intervals_neoZ$Chr)


p <- ggplot(Mpv2.full.map,
       aes(x = V4/1000000, 
           y = -log10(Mpv2.full.q.values),color=pcat)) + 
  geom_point(size = 0.4,shape = 19)+
  geom_segment(data=intervals_df,aes(x=start/1e6,xend=end/1e6),y=y_top,yend=y_top,color='blue',linewidth=2,inherit.aes=FALSE)+
  geom_segment(data = intervals_neoZ,aes(x=start,xend=end),y = y_top, yend = y_top,color = 'black', linewidth = 2.8,inherit.aes=FALSE) +
  geom_segment(data=intervals_neoZ,aes(x=start,xend=end),y=y_top,yend=y_top,color=c('black','white','#857E13'),linewidth=2,inherit.aes=FALSE)+
  facet_wrap(.~Chr,scales="free_x")+
  geom_hline(yintercept = -log10(0.005),linetype = "dashed",color = "darkred") +
  theme_minimal() +
  labs(
    x = "Chromosome position (Mb)",
    y = expression(-log[10](italic(p)))
  )+
  scale_color_manual(values=c("grey","darkred"))+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="none",
        panel.grid.major = element_line(color = "grey99", linewidth = 0.5),
        panel.grid.minor = element_line(color = "white", linewidth = 0.25),
        panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.5))




p


# Create dummy legend
######################################################################################
# Create a tiny data frame with one row per legend category
# Example data
# Ensure segment_type is a factor
df_legend <- data.frame(
  x = 1:4,      
  y = 1:4,      
  segment_type = factor(c("ancestral Z", "added Z", "new PAR", "MDS outlier regions"),
                        levels = c("ancestral Z", "added Z", "new PAR", "MDS outlier regions"))
)

# Generate the dummy plot to extract the legend
dummy_legend_plot <- ggplot(df_legend, aes(x, y, fill = segment_type, color = segment_type)) +
  geom_tile(width = 0.7, height = 0.3, size = 1) +  # Rectangular legend keys
  scale_fill_manual(
    name = "Segments",  # Ensure a legend title is set
    values = c(
      "ancestral Z" = "black",
      "added Z"     = "white",
      "new PAR"     = "#857E13",
      "MDS outlier regions"        = "blue"
    )
  ) +
  scale_color_manual(  
    name = "Segments",  # Match the name to avoid splitting into two legends
    values = c(
      "ancestral Z" = "black",
      "added Z"     = "black",  # Black outline for white fill
      "new PAR"     = "#857E13",
      "MDS outlier regions"        = "blue"
    )
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(
        color = c("black", "black", "#857E13", "blue"),  # Black outline for white tile
        fill  = c("black", "white", "#857E13", "blue"),  # Correct interior colors
        size  = 3  # Adjust size for better visibility
      )
    )
  ) +
  theme_void() +
  theme(legend.position = "top", legend.title = element_blank())

# Check again if the legend exists
print(dummy_legend_plot)

components <- ggplotGrob(dummy_legend_plot)

# Print all grobs in the plot to see where the legend is
print(components)
custom_legend <- gtable::gtable_filter(components, "guide-box-top")

# Print the extracted legend (This is now a gtable, NOT a full ggplot)
grid::grid.draw(custom_legend)



final_plot <- plot_grid(
  custom_legend,   # Extracted legend
  p,               # Your main ggplot figure
  ncol = 1,        # Arrange them vertically
  rel_heights = c(0.1,1)  # Give more space to the main plot
)

print(final_plot)





# Save the plot as a PNG with 400 dpi, for example 10 x 7 inches:
ggsave("Mpv2.full.map_PC1_2_legend_updated.png", plot = final_plot, dpi = 350, width = 14, height = 10, units = "in")

################################################################################
# Create a list of all SNPs and their locations with adjusted p-value <0.005 (to be moderately conservative)


# combine map with adjusted p-values
Mpv2.full.map <- cbind(Mpv2.full.q.values,Mmap)

MsigSNPs <- Mpv2.full.map[Mpv2.full.map$Mpv2.full.q.values<0.005,]


write.table(MsigSNPs,file="Mpv2.full.map.sigSNPs_0.005_threshold.txt",quote=FALSE,row.names=FALSE,sep="\t")



################################################################################
# Plot just the neo-Z with gene annotations of significant hits


#write.table(neoZ_results,"neoZ_results_GEA.txt",quote=FALSE,row.names = FALSE)
# added gene names manually to these results

#write.table(neoZ_results,"neoZ_results_GEA.txt",quote=FALSE,row.names = FALSE)
# added gene names manually to these results

GEA_neoZ <- read.table('neoZ_results_GEA.txt',header=T)
GEA_neoZ_filt <- GEA_neoZ[GEA_neoZ$gene!='unknown',]

GEA_neoZ_filt_order <- GEA_neoZ_filt[order(GEA_neoZ_filt$gene),]
GEA_neoZ_filt_order_nondup <- GEA_neoZ_filt_order[!duplicated(GEA_neoZ_filt_order$gene),]

Mpv2.neoZ.map <- Mpv2.full.map[Mpv2.full.map$V1=="scaffold_2_RagTag",]

# position variable for geom_segment()
y_top <- max(-log10(Mpv2.neoZ.map$Mpv2.full.q.values), na.rm = TRUE) * 1.1

p <- ggplot(Mpv2.neoZ.map,
            aes(x = V4/1000000, 
                y = -log10(Mpv2.full.q.values),color=pcat)) + 
  geom_point(size = 0.4,shape = 19)+
  geom_segment(data = intervals_neoZ,aes(x=start,xend=end),y = y_top, yend = y_top,color = 'black', linewidth = 2.8,inherit.aes=FALSE) +
  geom_segment(data=intervals_neoZ,aes(x=start,xend=end),y=y_top,yend=y_top,color=c('black','white','#857E13'),linewidth=2,inherit.aes=FALSE)+
  geom_hline(yintercept = -log10(0.005),linetype = "dashed",color = "darkred") +
  theme_minimal() +
  ylim(0,7)+
  labs(
    x = "Chromosome position (Mb)",
    y = expression(-log[10](italic(p)))
  )+
  scale_color_manual(values=c("grey","darkred"))+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="none",
        panel.grid.major = element_line(color = "grey99", linewidth = 0.5),
        panel.grid.minor = element_line(color = "white", linewidth = 0.25),
        panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.5))+
  geom_text_repel(
    data = GEA_neoZ_filt_order_nondup,       # Separate annotation data frame
    aes(x = V4/1000000, y = -log10(Mpv2.full.q.values), label = gene),
    color = "blue",        # Customize text color
    vjust = -1,            # Adjust vertical position
    size = 5,
    max.overlaps = 20
  )





p

ggsave("Mpv2.neoZ.gene_labels_significant_hits.png", plot = p, dpi = 350, width = 10, height = 6, units = "in")



################################################################################
# full=TRUE, K=3
Mpv3.full.q.values <- p.adjust(Mpv3.full$pvalues,method="BH")
plot(-log10(Mpv3.full.q.values), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.00005), lty = 2, col = "darkred")
points(ChrZ_target, -log10(Mpv3.full.q.values[ChrZ_target]), col = "darkred")
points(Chr4_target, -log10(Mpv3.full.q.values[Chr4_target]), col = "red")
points(Chr6_target, -log10(Mpv3.full.q.values[Chr6_target]), col = "blue")
points(Chr8_target, -log10(Mpv3.full.q.values[Chr8_target]), col = "lightblue")
points(Chr9_target, -log10(Mpv3.full.q.values[Chr9_target]), col = "magenta")
points(Chr10_target, -log10(Mpv3.full.q.values[Chr10_target]), col = "orange")
points(Chr13_target, -log10(Mpv3.full.q.values[Chr13_target]), col = "green")
points(Chr14_target, -log10(Mpv3.full.q.values[Chr14_target]), col = "yellow")
points(Chr15_target, -log10(Mpv3.full.q.values[Chr15_target]), col = "black")
points(Chr16_target, -log10(Mpv3.full.q.values[Chr16_target]), col = "darkgreen")
points(Chr17_target, -log10(Mpv3.full.q.values[Chr17_target]), col = "pink")
points(Chr18_target, -log10(Mpv3.full.q.values[Chr18_target]), col = "brown")
points(Con15l_target, -log10(Mpv3.full.q.values[Con15l_target]), col = "darkblue")
points(Con45l_target, -log10(Mpv3.full.q.values[Con45l_target]), col = "darkorange")

# Add a legend with filled points
legend("topleft", 
       legend = c("ChrZ", "Chr4", "Chr6", "Chr8", "Chr9", "Chr10", 
                  "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", 
                  "Chr18", "Con15l", "Con45l"),
       col = c("darkred", "red", "blue", "lightblue", "magenta", "orange", 
               "green", "yellow", "black", "darkgreen", "pink", 
               "brown", "darkblue", "darkorange"),
       pch = 16,  # Use filled circles
       bty = "n") # Remove box around legend






########################################################################
#K=2
Mpv2.PC1.q.values <- p.adjust(Mpv2$pvalues[1,],method="BH")
plot(-log10(Mpv2.PC1.q.values), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.00005), lty = 2, col = "darkred")
points(ChrZ_target, -log10(Mpv2.PC1.q.values[ChrZ_target]), col = "darkred")

Mpv2.PC2.q.values <- p.adjust(Mpv2$pvalues[2,],method="BH")
plot(-log10(Mpv2.PC2.q.values), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.00005), lty = 2, col = "darkred")
points(ChrZ_target, -log10(Mpv2.PC2.q.values[ChrZ_target]), col = "darkred")



#K=3
Mpv3.PC1.q.values <- p.adjust(Mpv3$pvalues[1,],method="BH")
plot(-log10(Mpv3.PC1.q.values), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.00005), lty = 2, col = "darkred")
points(ChrZ_target, -log10(Mpv3.PC1.q.values[ChrZ_target]), col = "darkred")

Mpv3.PC2.q.values <- p.adjust(Mpv3$pvalues[2,],method="BH")
plot(-log10(Mpv3.PC2.q.values), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.00005), lty = 2, col = "darkred")
points(ChrZ_target, -log10(Mpv3.PC2.q.values[ChrZ_target]), col = "darkred")


###########################################################################
# Combine the neo-Z and the autosomes into one plot and explain in methods
###########################################################################

# combine map with adjusted p-values
Mpv2.full.map <- cbind(Mpv2.full.q.values,Mmap)
pv2.full.map <- cbind(pv2.full.q.values,map)


Mpv2.neoZ <- Mpv2.full.map[Mpv2.full.map$V1=="scaffold_2_RagTag",]
colnames(Mpv2.neoZ)[1] <- 'q.values'
colnames(pv2.full.map)[1] <- 'q.values'

combined <- rbind(pv2.full.map,Mpv2.neoZ)


# Updated SV breakpoints 26 Feb 2025
lines_df <- data.frame(
  Chr = c(
    # "Chr_4", "Chr_4",
    "Chr_4", "Chr_4",
    "Chr_6", "Chr_6",
    # "Chr_6", "Chr_6",
    "Chr_8", "Chr_8",
    "Chr_8", "Chr_8",
    # "Chr_9", "Chr_9",
    "Chr_9", "Chr_9",
    "pri#ptg000015l", "pri#ptg000015l",
    "Chr_10", "Chr_10",
    "Chr_13", "Chr_13",
    "Chr_14", "Chr_14",
    "Chr_15", "Chr_15",
    "Chr_17", "Chr_17",
    "Chr_18", "Chr_18",
    "pri#ptg000045l", "pri#ptg000045l"
  ),
  vline_pos = c(
    # Chr_4_:95375035-105124963
    # 95375035, 105124963,
    # Chr_4_:107824981-120434650
    # 107824981, 120434650,
    94877058, 120434650, 
    # Chr_6_:8125001-14575038
    #8125001, 14575038,
    # Chr_6_:15924948-328749905
    #15924948, 32874990,
    8125001,32874990,
    # Chr_8_:1-10724985
    1, 10724985,
    # Chr_8_:5597491-37874992
    56196169, 37874992,
    # Chr_9_:16374481-25249745
    #16374481, 2524974,
    # Chr_9_:16374481-18339159
    1631342, 18339159,
    # pri#ptg000015l:6864979-17909547
    6864979, 17909547,
    # Chr_10_:9675047-22345924
    9675047, 22345924,
    # Chr_13_:15724920-25076809
    # Chr_13_:14500000-25076809
    15000000, 23500000,
    # Chr_14_:1-101749905
    1, 10174990,
    # Chr_15_:1-877497
    1, 8774970,
    # Chr_17_:12914902-19617658
    12708285, 19617658,
    # Chr_18_:1865069-10487981
    2865069, 10487981,
    # pri#ptg000045l:1-128348605
    1, 13100247
  ),
  stringsAsFactors = TRUE
)


combined$Chr <- as.character(combined$V1)

lines_df$Chr <- as.character(lines_df$Chr)

intervals_df <- lines_df %>%
  group_by(Chr) %>%
  arrange(vline_pos) %>%
  # Create pairs: rows 1-2 form the first interval, rows 3-4 form the second, etc.
  mutate(row_index = row_number()) %>%
  # For each pair, the "start" is an odd row_index, the "end" is the following even row_index
  # We'll pivot those out into one row per pair.
  summarize(
    start = vline_pos[row_index %% 2 == 1], 
    end   = vline_pos[row_index %% 2 == 0],
    .groups = "keep"
  ) %>%
  # It's often handy to reorder columns or rename them
  select(Chr, start, end)



combined$in_interval <- FALSE

combined$Chr <- str_replace(combined$Chr,'scaffold','Chr')
combined$Chr <- str_remove(combined$Chr,'_RagTag')

# For each interval, mark those points as in_interval
for(i in seq_len(nrow(intervals_df))) {
  chr_i   <- intervals_df$Chr[i]
  start_i <- intervals_df$start[i]
  end_i   <- intervals_df$end[i]
  
  combined$in_interval <- combined$in_interval | 
    (combined$Chr == chr_i & combined$V4 >= start_i & combined$V4 <= end_i)
}



table(combined$Chr)

combined$Chr[combined$V1=="scaffold_2_RagTag"]<-"Neo_Z"

combined$Chr <- factor(combined$Chr,levels=c("Neo_Z","Chr_1","Chr_3","Chr_4","Chr_6",
                                                       "Chr_7","Chr_8","Chr_9","Chr_10",
                                                       "Chr_11","Chr_12","Chr_13","Chr_14",
                                                       "Chr_15","Chr_16","Chr_17","Chr_18",
                                                       "Chr_19","Chr_20","Chr_21","Chr_22",
                                                       "Chr_23","Chr_24","Chr_25","Chr_26",
                                                       "Chr_27","Chr_28","Chr_29","Chr_30",
                                                       "Chr_31","Chr_32","Chr_33","Chr_34",
                                                       "Chr_35","Chr_36","Chr_37","Chr_38",
                                                       "Chr_39","pri#ptg000015l","pri#ptg000031l",
                                                       "pri#ptg000045l","pri#ptg000050l","pri#ptg000051l",
                                                       "pri#ptg000085l","pri#ptg000090l"))




lines_df$Chr <- factor(lines_df$Chr,levels=c("Neo_Z","Chr_1","Chr_3","Chr_4","Chr_6",
                                             "Chr_7","Chr_8","Chr_9","Chr_10",
                                             "Chr_11","Chr_12","Chr_13","Chr_14",
                                             "Chr_15","Chr_16","Chr_17","Chr_18",
                                             "Chr_19","Chr_20","Chr_21","Chr_22",
                                             "Chr_23","Chr_24","Chr_25","Chr_26",
                                             "Chr_27","Chr_28","Chr_29","Chr_30",
                                             "Chr_31","Chr_32","Chr_33","Chr_34",
                                             "Chr_35","Chr_36","Chr_37","Chr_38",
                                             "Chr_39","pri#ptg000015l","pri#ptg000031l",
                                             "pri#ptg000045l","pri#ptg000050l","pri#ptg000051l",
                                             "pri#ptg000085l","pri#ptg000090l"))


intervals_df$Chr <- factor(intervals_df$Chr,levels=c("Neo_Z","Chr_1","Chr_3","Chr_4","Chr_6",
                                                     "Chr_7","Chr_8","Chr_9","Chr_10",
                                                     "Chr_11","Chr_12","Chr_13","Chr_14",
                                                     "Chr_15","Chr_16","Chr_17","Chr_18",
                                                     "Chr_19","Chr_20","Chr_21","Chr_22",
                                                     "Chr_23","Chr_24","Chr_25","Chr_26",
                                                     "Chr_27","Chr_28","Chr_29","Chr_30",
                                                     "Chr_31","Chr_32","Chr_33","Chr_34",
                                                     "Chr_35","Chr_36","Chr_37","Chr_38",
                                                     "Chr_39","pri#ptg000015l","pri#ptg000031l",
                                                     "pri#ptg000045l","pri#ptg000050l","pri#ptg000051l",
                                                     "pri#ptg000085l","pri#ptg000090l"))




combined$Chr
# plot

# Create a categorical variable for SNPs that are significant
combined$pcat <- FALSE
combined$pcat[combined$q.values<0.005] <- TRUE


# position variable for geom_segment()
y_top <- max(-log10(combined$q.values), na.rm = TRUE) * 1.01

# Plot intervals for neo-Z
intervals_neoZ <- data.frame(Chr=c("Neo_Z","Neo_Z","Neo_Z"),start=c(0,74.4,111.6),end=c(74.4,111.6,128.9))
intervals_neoZ$Chr <- as.factor(intervals_neoZ$Chr)


p <- ggplot(combined,
            aes(x = V4/1000000, 
                y = -log10(q.values),color=pcat)) + 
  geom_point(size = 0.4,shape = 19)+
  geom_segment(data=intervals_df,aes(x=start/1e6,xend=end/1e6),y=y_top,yend=y_top,color='blue',linewidth=2,inherit.aes=FALSE)+
  # Black outline segments (slightly thicker)
  geom_segment(data = intervals_neoZ,aes(x=start,xend=end),y = y_top, yend = y_top,color = 'black', linewidth = 2.8,inherit.aes=FALSE) +
  geom_segment(data=intervals_neoZ,aes(x=start,xend=end),y=y_top,yend=y_top,color=c('black','white','#857E13'),linewidth=2,inherit.aes=FALSE)+
  facet_wrap(.~Chr,scales="free_x")+
  geom_hline(yintercept = -log10(0.005),linetype = "dashed",color = "darkred") +
  theme_minimal() +
  labs(
    x = "Chromosome position (Mb)",
    y = expression(-log[10](italic(p)))
  )+
  scale_color_manual(values=c("grey","darkred"))+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="none",
        panel.grid.major = element_line(color = "grey99", linewidth = 0.5),
        panel.grid.minor = element_line(color = "white", linewidth = 0.25),
        panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.5))

p


################################################################
# just the inversion
inversion <- combined[combined$V1=="pri#ptg000015l",]
inversion$Chr <- as.character(inversion$Chr)
# position variable for geom_segment()
y_top <- max(-log10(inversion$q.values), na.rm = TRUE) * 1.01

p <- ggplot(inversion,
            aes(x = V4/1000000, 
                y = -log10(q.values),color=pcat)) + 
  geom_point(size = 1,shape = 19)+
 geom_segment(x=6864979/1e6, xend=17909547/1e6,y=y_top,yend=y_top,color='blue',linewidth=2)+
  # Black outline segments (slightly thicker)
  geom_hline(yintercept = -log10(0.005),linetype = "dashed",color = "darkred") +
  theme_minimal() +
  labs(
    x = "Chromosome position (Mb)",
    y = expression(-log[10](italic(p)))
  )+
  scale_color_manual(values=c("grey","darkred"))+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="none",
        panel.grid.major = element_line(color = "grey99", linewidth = 0.5),
        panel.grid.minor = element_line(color = "white", linewidth = 0.25),
        panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.5))

p





################################################################
# just the neo-Z


neoZ <- combined[combined$V1=="scaffold_2_RagTag",]
neoZ$Chr <- as.character(neoZ$Chr)
# position variable for geom_segment()
y_top <- max(-log10(neoZ$q.values), na.rm = TRUE) * 1.05

p <- ggplot(neoZ,
            aes(x = V4/1000000, 
                y = -log10(q.values),color=pcat)) + 
  geom_point(size = 1,shape = 19)+
  # Black outline segments (slightly thicker)
  geom_segment(data = intervals_neoZ,aes(x=start,xend=end),y = y_top, yend = y_top,color = 'black', linewidth = 2.8,inherit.aes=FALSE) +
  geom_segment(data=intervals_neoZ,aes(x=start,xend=end),y=y_top,yend=y_top,color=c('black','white','#857E13'),linewidth=2,inherit.aes=FALSE)+
  # Black outline segments (slightly thicker)
  geom_hline(yintercept = -log10(0.005),linetype = "dashed",color = "darkred") +
  theme_minimal() +
  labs(
    x = "Chromosome position (Mb)",
    y = expression(-log[10](italic(p)))
  )+
  ylim(0,7)+
  scale_color_manual(values=c("grey","darkred"))+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="none",
        panel.grid.major = element_line(color = "grey99", linewidth = 0.5),
        panel.grid.minor = element_line(color = "white", linewidth = 0.25),
        panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.5))

p









# Create dummy legend
######################################################################################
# Create a tiny data frame with one row per legend category
# Example data
# Ensure segment_type is a factor
df_legend <- data.frame(
  x = 1:4,      
  y = 1:4,      
  segment_type = factor(c("ancestral Z", "added Z", "new PAR", "MDS outlier regions"),
                        levels = c("ancestral Z", "added Z", "new PAR", "MDS outlier regions"))
)

# Generate the dummy plot to extract the legend
dummy_legend_plot <- ggplot(df_legend, aes(x, y, fill = segment_type, color = segment_type)) +
  geom_tile(width = 0.7, height = 0.3, size = 1) +  # Rectangular legend keys
  scale_fill_manual(
    name = "Segments",  # Ensure a legend title is set
    values = c(
      "ancestral Z" = "black",
      "added Z"     = "white",
      "new PAR"     = "#857E13",
      "MDS outlier regions"        = "blue"
    )
  ) +
  scale_color_manual(  
    name = "Segments",  # Match the name to avoid splitting into two legends
    values = c(
      "ancestral Z" = "black",
      "added Z"     = "black",  # Black outline for white fill
      "new PAR"     = "#857E13",
      "MDS outlier regions"        = "blue"
    )
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(
        color = c("black", "black", "#857E13", "blue"),  # Black outline for white tile
        fill  = c("black", "white", "#857E13", "blue"),  # Correct interior colors
        size  = 3  # Adjust size for better visibility
      )
    )
  ) +
  theme_void() +
  theme(legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(size=12))

# Check again if the legend exists
print(dummy_legend_plot)

components <- ggplotGrob(dummy_legend_plot)

# Print all grobs in the plot to see where the legend is
print(components)
custom_legend <- gtable::gtable_filter(components, "guide-box-top")

# Print the extracted legend (This is now a gtable, NOT a full ggplot)
grid::grid.draw(custom_legend)



final_plot <- plot_grid(
  custom_legend,   # Extracted legend
  p,               # Your main ggplot figure
  ncol = 1,        # Arrange them vertically
  rel_heights = c(0.1,1)  # Give more space to the main plot
)

print(final_plot)


ggsave("GEA_K2_PC1_PC2_all_autosomes_neoZ_breakpoints_legend_updated.png", 
       plot = final_plot, dpi = 350, width = 14, height = 10, units = "in")



