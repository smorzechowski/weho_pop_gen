# Run GEA analysis on White-eared Honeyeater resequencing data
# Sophia MacRae Orzechowski
# Use dummy variables for inversions
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
#ped2geno("Nleu_autos_lea_depth_filt_imputed_thin_sans_inversions_with_dummy_midpoint_plink.ped")
#ped2geno("Nleu_autos_neoZ_lea_males_depth_filt_imputed_thin_sans_inversions_with_dummy_midpoint_plink.ped")
#ped2geno("Nleu_autos_neoZ_lea_males_depth_filt_imputed_thin_sans_inversions_neoZ_with_dummy_midpoint_plink.ped")


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

################################################################################
## All autosomes except new PAR region, 65 individuals


#gen.imp <- read.geno("Nleu_autos_lea_depth_filt_imputed_thin_sans_inversions_with_dummy_plink.geno")
gen.imp <- read.geno("Nleu_autos_lea_depth_filt_imputed_thin_sans_inversions_with_dummy_midpoint_plink.geno")
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


dim(pred)
str(pred)
class(pred) 



#map <- read.table("Nleu_autos_lea_depth_filt_imputed_thin_sans_inversions_with_dummy_plink.map",header=F,comment.char = "")
map <- read.table("Nleu_autos_lea_depth_filt_imputed_thin_sans_inversions_with_dummy_midpoint_plink.map",header=F,comment.char = "")

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
pv2.full <- lfmm2.test(object = mod2,
                 input = gen.imp,
                 env = pred,
                 full=TRUE)

pv3.full <- lfmm2.test(object = mod3,
                 input = gen.imp,
                 env = pred,
                 full=TRUE)

# Run univariate tests for each variable, NOT accounting for shared variance (diagonal covariance matrix)
pv2 <- lfmm2.test(object = mod2,
                       input = gen.imp,
                       env = pred,
                       full=FALSE)

pv3 <- lfmm2.test(object = mod3,
                       input = gen.imp,
                       env = pred,
                       full=FALSE)

pv4 <- lfmm2.test(object = mod4,
                  input = gen.imp,
                  env = pred,
                  full=FALSE)

################################################################################
##### Built-in FDR procedure #####


# Calculate adjusted pvalues for the full model
pv2.full.q.values <- p.adjust(pv2.full$pvalues, method = "BH")
pv3.full.q.values <- p.adjust(pv3.full$pvalues, method = "BH")


################################################################################
#Create plots
################################################################################
# Full: all variables combined for both K=2 and K=3
# ggplot2 version: Plot of all variables combined with BH correction K=2

#pv2.full.map <- cbind(pv2.full.q.values,map)
pv2.full.map.dummy <- cbind(pv2.full.q.values,map)
hits <- data.frame(table(pv2.full.map.dummy$V1[pv2.full.map.dummy$pv2.full.q.values<0.005]))


source("/n/netscratch/edwards_lab/Lab/smorzechowski/meliphagid/analysis/2024-11-03/19-lea/chromosome_limits.R")

dummy_limits <- dummy_limits %>% 
  mutate(x_dummy = x_dummy / 1e6)


# plot

pv2.full.map$CHROM <- pv2.full.map$V1

ggplot(pv2.full.map,
       aes(x = V4/1000000, 
           y = -log10(pv2.full.q.values))) + 
  # 1) Plot ALL points in grey
  geom_point(size = 0.4,# roughly corresponds to cex = 0.4 in base R
             shape = 19)+
  geom_blank(data = dummy_limits, aes(x = x_dummy, y = y_dummy)) +
  facet_wrap(.~CHROM,scales="free_x")+
  geom_hline(
    yintercept = -log10(0.00005), 
    linetype = "dashed", 
    color = "darkred"
  )


# Create data frame for vertical lines 
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
    95375035, 120434650, 
    # Chr_6_:8125001-14575038
    #8125001, 14575038,
    # Chr_6_:15924948-328749905
    #15924948, 32874990,
    8125001,32874990,
    # Chr_8_:1-10724985
    1, 10724985,
    # Chr_8_:5597491-37874992
    55974910, 37874992,
    # Chr_9_:16374481-25249745
    #16374481, 2524974,
    # Chr_9_:16374481 updated -18339159
    1637448, 18339159,
    # pri#ptg000015l:6864979-17909547
    6864979, 17909547,
    # Chr_10_:9675047-22345924
    9675047, 22345924,
    # Chr_13_:15724920-25076809
    # Chr_13_:14500000-25076809
    15000000, 25076809,
    # Chr_14_:1-101749905
    1, 10174990,
    # Chr_15_:1-877497
    1, 8774970,
    # Chr_17_:12914902-19617658
    12914902, 19617658,
    # Chr_18_:1865069-10487981
    2865069, 10487981,
    # pri#ptg000045l:1-128348605
    1, 12834860
  ),
  stringsAsFactors = FALSE
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



pv2.full.map.dummy$in_interval <- FALSE

pv2.full.map.dummy$Chr <- str_replace(pv2.full.map.dummy$V1,'scaffold','Chr')
pv2.full.map.dummy$Chr <- str_remove(pv2.full.map.dummy$Chr,'_RagTag')

# For each interval, mark those points as in_interval
for(i in seq_len(nrow(intervals_df))) {
  chr_i   <- intervals_df$Chr[i]
  start_i <- intervals_df$start[i]
  end_i   <- intervals_df$end[i]
  
  pv2.full.map.dummy$in_interval <- pv2.full.map.dummy$in_interval | 
    (pv2.full.map.dummy$Chr == chr_i & pv2.full.map.dummy$V4 >= start_i & pv2.full.map.dummy$V4 <= end_i)
}


table(pv2.full.map.dummy$Chr)

pv2.full.map.dummy$Chr <- factor(pv2.full.map.dummy$Chr,levels=c("Chr_1","Chr_3","Chr_4","Chr_6",
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



#pv2.full.map.dummy$Chr <- as.character(pv2.full.map.dummy$Chr)
#intervals_df$Chr <- as.character(intervals_df$Chr)

str(pv2.full.map.dummy$Chr)
table(pv2.full.map.dummy$Chr)
# plot

dummy_limits$Chr <- dummy_limits$CHROM
dummy_limits$Chr <- str_replace(dummy_limits$Chr,'scaffold','Chr')
dummy_limits$Chr <- str_remove(dummy_limits$Chr,'_RagTag')

dummy_limits$Chr <- factor(dummy_limits$Chr,levels=c("Chr_1","Chr_3","Chr_4","Chr_6",
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


lines_df$Chr <- factor(lines_df$Chr,levels=c("Chr_1","Chr_3","Chr_4","Chr_6",
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


intervals_df$Chr <- factor(intervals_df$Chr,levels=c("Chr_1","Chr_3","Chr_4","Chr_6",
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




# Create dummy legend
######################################################################################
# Create a tiny data frame with one row per legend category
# Example data
# Ensure segment_type is a factor

# Create the dummy legend data frame
df_legend <- data.frame(
  x = c(1, 1),      
  y = c(1, 2),      
  segment_type = factor(c("MDS outlier regions", "dummy SNP"), levels = c("MDS outlier regions", "dummy SNP"))
)

# Generate the dummy plot with both segments & the SNP circle
dummy_legend_plot <- ggplot(df_legend, aes(x, y, fill = segment_type, color = segment_type, shape = segment_type)) +
  # Rectangular legend item for "PCRs"
  geom_tile(data = df_legend[df_legend$segment_type == "MDS outlier regions", ], width = 0.7, height = 0.3) +
  
  # Large red circle legend item for "dummy SNP"
  geom_point(data = df_legend[df_legend$segment_type == "dummy SNP", ], 
             size = 6, stroke = 1.5) +
  
  # Define manual scales for color, fill, and shape
  scale_fill_manual(
    name = "Legend",  # Legend title
    values = c("MDS outlier regions" = "blue", "dummy SNP" = "red")  # Both must be defined here
  ) +
  scale_color_manual(
    name = "Legend",
    values = c("MDS outlier regions" = "blue", "dummy SNP" = "red")  # Colors must be consistent
  ) +
  scale_shape_manual(
    name = "Legend",
    values = c("MDS outlier regions" = NA, "dummy SNP" = 1)  # Shape 21 is a circle
  ) +
  
  # Force ggplot to combine them into a single legend
  guides(
    fill = guide_legend(override.aes = list(shape = c(NA, 1), color = c("blue", "red"), size = c(NA, 6))),
    color = "none",  # Hide redundant color legend
    shape = "none"   # Hide redundant shape legend
  ) +
  
  theme_void() +
  theme(legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(size=12))

# Display the legend
print(dummy_legend_plot)

# Extract the legend
components <- ggplotGrob(dummy_legend_plot)

# Print all grobs in the plot to see where the legend is
print(components)
custom_legend <- gtable::gtable_filter(components, "guide-box")

# Print the extracted legend (This is now a gtable, NOT a full ggplot)
grid::grid.draw(custom_legend)

#################################################################################
# Create a categorical variable for SNPs that are significant
pv2.full.map.dummy$pcat <- FALSE
pv2.full.map.dummy$pcat[pv2.full.map.dummy$pv2.full.q.values<0.005] <- TRUE

table(pv2.full.map.dummy$pcat)



y_top <- max(-log10(pv2.full.map.dummy$pv2.full.q.values), na.rm = TRUE) * 1.1




p <- ggplot(pv2.full.map.dummy,
       aes(x = V4/1000000, 
           y = -log10(pv2.full.q.values),color=pcat,size=in_interval)) + 
  # 1) Plot ALL points in grey
  geom_point(size = 0.4,# roughly corresponds to cex = 0.4 in base R
    shape = 19
  ) +
#geom_vline(
#    data = lines_df,
#    aes(xintercept = vline_pos / 1e6),
#    color = "blue",
#    linetype = "dashed",
#  )+
  geom_segment(data=intervals_df,aes(x=start/1e6,xend=end/1e6),y=y_top,yend=y_top,color='blue',linewidth=2,inherit.aes=FALSE)+
  geom_blank(data = dummy_limits, aes(x = x_dummy, y = y_dummy),inherit.aes=FALSE) +
  facet_wrap(.~Chr,scales="free_x")+
  ylim(0,12)+
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
 # scale_color_manual(values=c("darkblue","darkblue"))+
  scale_color_manual(values=c("grey","darkred"))+
  theme(legend.position = "none",
        axis.title=element_text(size=15),
        strip.text = element_text(size=13),
        axis.text= element_text(size=12),
        panel.grid.major = element_line(color = "grey99", linewidth = 0.5),
        panel.grid.minor = element_line(color = "white", linewidth = 0.25),
        panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.5))+
  scale_size_manual(values=c(1,15))+
  geom_point(data = subset(pv2.full.map.dummy, in_interval == TRUE),
             shape = 1,       # shape 1 is an unfilled circle
             size = 4,        # adjust size as needed; this is the size of the outline circle
             stroke = 1,      # thickness of the circle's outline
             color = "red")   # color of the circle outline



p



final_plot <- plot_grid(
  custom_legend,   # Extracted legend
  p,               # Your main ggplot figure
  ncol = 1,        # Arrange them vertically
  rel_heights = c(0.1,1)  # Give more space to the main plot
)

print(final_plot)


ggsave("GEA_K2_PC1_PC2_all_autosomes_dummy_PCR_breakpoints_legend_updated.png", 
       plot = final_plot, dpi = 350, width = 14, height = 10, units = "in")





################################################################################
## All autosomes and neo-Z, including new PAR, 43 males
################################################################################



#Mgen.imp <- read.geno("Nleu_autos_neoZ_lea_males_depth_filt_imputed_thin_sans_inversions_with_dummy_midpoint_plink.geno")
Mgen.imp <- read.geno("Nleu_autos_neoZ_lea_males_depth_filt_imputed_thin_sans_inversions_neoZ_with_dummy_midpoint_plink.geno")
dim(Mgen.imp)


#Mpred <- read.env("Nleu_bioclim_variables_43males_autos_neoZ_ordered_noheader.env")
Mpred <- read.env("Nleu_bioclim_PC1_PC2_43males_autos_neoZ.env")

#Mpred <- Mpred[,c(19,8)]

dim(Mgen.imp)
dim(Mpred)
str(Mpred)
class(Mpred) 

#Mmap <- read.table("Nleu_autos_neoZ_lea_males_depth_filt_imputed_thin_sans_inversions_with_dummy_midpoint_plink.map",header=F,comment.char = "")
Mmap <- read.table("Nleu_autos_neoZ_lea_males_depth_filt_imputed_thin_sans_inversions_neoZ_with_dummy_midpoint_plink.map",header=F,comment.char = "")
Mmod2 <- lfmm2(input = Mgen.imp, env = Mpred, K = 2)


Mmod3 <- lfmm2(input = Mgen.imp, env = Mpred, K = 3)

################################################################################
# evaluate collinearity between latent factors and environmental variables
# Extract the factor matrix U (n x K)

U <- Mmod2@U 

# Evaluate collinearity for each environment variable:
for (j in 1:ncol(Mpred)) {
  lm_fit <- lm(Mpred[,j] ~ U)
  r2_val <- summary(lm_fit)$r.squared
  cat(colnames(Mpred)[j], "-> R^2 with factors:", round(r2_val, 3), "\n")
}


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


#######################################################################
# Calculate adjusted p-values
# full = TRUE, K=2
Mpv2.full.q.values <- p.adjust(Mpv2.full$pvalues,method="BH")
Mpv3.full.q.values <- p.adjust(Mpv3.full$pvalues,method="BH")


######################################################################
#Create plots
######################################################################
# ggplot2 version: Plot of all variables combined with BH correction K=2


Mpv2.full.map <- cbind(Mpv2.full.q.values,Mmap)

Mpv2.full.map <- cbind(Mpv3.full.q.values,Mmap)

source("/n/netscratch/edwards_lab/Lab/smorzechowski/meliphagid/analysis/2024-11-03/19-lea/chromosome_limits.R")

dummy_limits <- dummy_limits %>% 
  mutate(x_dummy = x_dummy / 1e6)


# plot

Mpv2.full.map$CHROM <- Mpv2.full.map$V1

ggplot(Mpv2.full.map,
       aes(x = V4/1000000, 
           y = -log10(Mpv2.full.q.values))) + 
  # 1) Plot ALL points in grey
  geom_point(size = 0.4,# roughly corresponds to cex = 0.4 in base R
             shape = 19)+
#  geom_blank(data = dummy_limits, aes(x = x_dummy, y = y_dummy)) +
  facet_wrap(.~CHROM,scales="free_x")
#  geom_hline(
#    yintercept = -log10(0.00005), 
#    linetype = "dashed", 
#    color = "darkred"
#  )


# Create data frame for vertical lines 
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
    95375035, 120434650, 
    # Chr_6_:8125001-14575038
    #8125001, 14575038,
    # Chr_6_:15924948-328749905
    #15924948, 32874990,
    8125001,32874990,
    # Chr_8_:1-10724985
    1, 10724985,
    # Chr_8_:5597491-37874992
    55974910, 37874992,
    # Chr_9_:16374481-25249745
    #16374481, 2524974,
    # Chr_9_:16374481 updated -18339159
    1637448, 18339159,
    # pri#ptg000015l:6864979-17909547
    6864979, 17909547,
    # Chr_10_:9675047-22345924
    9675047, 22345924,
    # Chr_13_:15724920-25076809
    # Chr_13_:14500000-25076809
    15000000, 25076809,
    # Chr_14_:1-101749905
    1, 10174990,
    # Chr_15_:1-877497
    1, 8774970,
    # Chr_17_:12914902-19617658
    12914902, 19617658,
    # Chr_18_:1865069-10487981
    2865069, 10487981,
    # pri#ptg000045l:1-128348605
    1, 12834860
  ),
  stringsAsFactors = FALSE
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


table(Mpv2.full.map$Chr)

Mpv2.full.map$Chr <- factor(Mpv2.full.map$Chr,levels=c("Chr_1","Chr_3","Chr_4","Chr_6",
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



#Mpv2.full.map$Chr <- as.character(Mpv2.full.map$Chr)
#intervals_df$Chr <- as.character(intervals_df$Chr)

str(Mpv2.full.map$Chr)
table(Mpv2.full.map$Chr)
# plot

dummy_limits$Chr <- dummy_limits$CHROM
dummy_limits$Chr <- str_replace(dummy_limits$Chr,'scaffold','Chr')
dummy_limits$Chr <- str_remove(dummy_limits$Chr,'_RagTag')

dummy_limits$Chr <- factor(dummy_limits$Chr,levels=c("Chr_1","Chr_3","Chr_4","Chr_6",
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


lines_df$Chr <- factor(lines_df$Chr,levels=c("Chr_1","Chr_3","Chr_4","Chr_6",
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


intervals_df$Chr <- factor(intervals_df$Chr,levels=c("Chr_1","Chr_3","Chr_4","Chr_6",
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


y_top <- max(-log10(Mpv2.full.map$Mpv2.full.q.values), na.rm = TRUE) * 1.1


ggplot(Mpv2.full.map,
       aes(x = V4/1000000, 
           y = -log10(Mpv2.full.q.values),color=in_interval,size=in_interval)) + 
  # 1) Plot ALL points in grey
  geom_point(size = 0.4,# roughly corresponds to cex = 0.4 in base R
             shape = 19
  ) +
  #geom_vline(
  #    data = lines_df,
  #    aes(xintercept = vline_pos / 1e6),
  #    color = "blue",
  #    linetype = "dashed",
  #  )+
  geom_segment(data=intervals_df,aes(x=start/1e6,xend=end/1e6),y=y_top,yend=y_top,color='blue',linewidth=2,inherit.aes=FALSE)+
  geom_blank(data = dummy_limits, aes(x = x_dummy, y = y_dummy),inherit.aes=FALSE) +
  facet_wrap(.~Chr,scales="free_x")+
  ylim(0,12)+
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
  scale_color_manual(values=c("darkblue","darkblue"))+
  theme(legend.position = "none",
        strip.text = element_text(size=13),
        axis.text= element_text(size=12))+
  scale_size_manual(values=c(1,15))+
  geom_point(data = subset(Mpv2.full.map, in_interval == TRUE),
             shape = 1,       # shape 1 is an unfilled circle
             size = 4,        # adjust size as needed; this is the size of the outline circle
             stroke = 1,      # thickness of the circle's outline
             color = "red")   # color of the circle outline








##########################################################################
##########################################################################
##########################################################################

# Create data frame for vertical lines 
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
    95375035, 120434650, 
    # Chr_6_:8125001-14575038
    #8125001, 14575038,
    # Chr_6_:15924948-328749905
    #15924948, 32874990,
    8125001,32874990,
    # Chr_8_:1-10724985
    1, 10724985,
    # Chr_8_:5597491-37874992
    55974910, 37874992,
    # Chr_9_:16374481-25249745
    #16374481, 2524974,
    # Chr_9_:16374481 updated-18339159
    1637448, 18339159,
    # pri#ptg000015l:6864979-17909547
    6864979, 17909547,
    # Chr_10_:9675047-22345924
    9675047, 22345924,
    # Chr_13_:15724920-25076809
    # Chr_13_:14500000-25076809
    15000000, 25076809,
    # Chr_14_:1-101749905
    1, 10174990,
    # Chr_15_:1-877497
    1, 8774970,
    # Chr_17_:12914902-19617658
    12914902, 19617658,
    # Chr_18_:1865069-10487981
    2865069, 10487981,
    # pri#ptg000045l:1-128348605
    1, 12834860
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


table(Mpv2.full.map$Chr)

Mpv2.full.map$Chr <- factor(Mpv2.full.map$Chr,levels=c("Chr_1","Chr_2","Chr_3","Chr_4","Chr_6",
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

Mpv2.full.map$pcat <- FALSE
Mpv2.full.map$pcat[Mpv2.full.map$Mpv2.full.q.values<0.0005] <- TRUE

ggplot(Mpv2.full.map,
       aes(x = V4/1000000, 
           y = -log10(Mpv2.full.q.values),color=pcat)) + 
  geom_point(size = 0.4,shape = 19)+
  facet_wrap(.~Chr,scales="free_x")+
  geom_hline(yintercept = -log10(0.00005),linetype = "dashed",color = "darkred") +
  theme_minimal() +
  labs(
    x = "Chromosome position (Mb)",
    y = expression(-log[10](italic(p)))
  )+
  scale_color_manual(values=c("grey","darkblue"))+
  theme(legend.position = "none",
        strip.text = element_text(size=13))






