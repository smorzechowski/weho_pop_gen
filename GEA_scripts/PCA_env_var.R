# Run PCA on environmental variables from bioclim
# Create env matrix of the top PC axes
# Sophie MacRae Orzechowski
# January 2025

library(LEA)
library(factoextra)
library(pheatmap)
library(caret)
library(dplyr)
library(tidyr)
library(ggplot2)

getwd()

setwd("/n/netscratch/edwards_lab/Lab/smorzechowski/meliphagid/analysis/2024-11-03/19-lea/")

env_data <- read.table('Nleu_bioclim_variables_65ind_autos_ordered.txt',header=T)
#env_data <- read.table('Nleu_bioclim_variables_43males_autos_neoZ_ordered.txt',header=T)


# no longer need rev() because I have ordered from bio1-19
colnames(env_data) <- c("Annual_Mean_T","Mean_Diurnal_Range","Isothermality","T_Seasonality","Max_T_Warmest_Mo","Min_T_Coldest_Mo",
                       "T_Annual_Range","Mean_T_Wettest_Q","Mean_T_Driest_Q","Mean_T_Warmest_Q","Mean_T_Coldest_Q","Annual_Precip",
                       "Precip_Wettest_Mo","Precip_Driest_Mo","Precip_Seasonality","Precip_Wettest_Q","Precip_Driest_Q","Precip_Warmest_Q","Precip_Coldest_Q")
summary(env_data)


pca_res <- prcomp(env_data, center = TRUE, scale. = TRUE)
summary(pca_res)
pca_res$rotation
head(pca_res$x)
plot(pca_res, type = "l",cex=1.4)

# Extract eigenvalues (variance explained)
eig_values <- pca_res$sdev^2  # Squared singular values give eigenvalues
eig_values <- eig_values[1:10]
eig_percent <- eig_values / sum(eig_values) * 100  # Convert to percentage


p <- fviz_eig(pca_res) +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=15))+
  xlab("Dimensions: PC axes")+
  geom_text(aes(x = seq_along(eig_percent), 
                y = eig_percent + 2,  # Adjust position slightly above bars
                label = sprintf("%.1f%%", eig_percent)), 
            size = 5, vjust = 0, fontface = "bold")

pc_scores <- pca_res$x[, 1:2]  # columns PC1, PC2

p

ggsave("Supplemental_fig_PCA_scree_plot.png", 
       plot = p, dpi = 300, width = 6, height = 6, units = "in")

#write.env(pc_scores,"Nleu_bioclim_PC1_PC2_65ind_autos.env")
#write.env(pc_scores,"Nleu_bioclim_PC1_PC2_43males_autos_neoZ.env")



###########################################################################
# The loading of each individual PC axis 


loading <- data.frame(pca_res$rotation)[,1:2]
loading$env <- factor(row.names(loading))

# Convert to long format
loading_long <- loading %>%
  pivot_longer(cols = starts_with("PC"), 
               names_to = "PC", 
               values_to = "Loading")

loading_long$env <- factor(loading_long$env,levels=row.names(loading))
p <- ggplot(loading_long,aes(env,Loading,fill=PC))+
  geom_bar(stat='identity',position='dodge')+
  coord_flip()+
  xlab('')+
  ylab('Loading value')+
  theme_minimal()+
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=15),
        legend.title=element_blank(),
        legend.text=element_text(size=13))

p
ggsave("Supplemental_fig_PCA_loading_values.png", 
       plot = p, dpi = 300, width = 8, height = 10, units = "in")


###########################################################################
# Correlation between the 19 bioclim variables

cor_matrix <- cor(env_data)


# Mask the upper triangle
cor_matrix[upper.tri(cor_matrix)] <- NA

# Reverse the order of the columns and rows
cor_matrix_reversed <- cor_matrix[rev(rownames(cor_matrix)), rev(colnames(cor_matrix))]

# Print the reversed correlation matrix
print(cor_matrix_reversed)



# Plot the heatmap
#pheatmap(cor_matrix, na_col = "black")  # White cells for NA


p <- pheatmap(
  cor_matrix_reversed,
  cluster_rows = FALSE,     # no row dendrogram
  cluster_cols = FALSE,     # no column dendrogram
  display_numbers = TRUE,   # show the correlation value in each cell
  number_format = "%.2f",   # optional: format numbers to 2 decimal places
  na_col="black",
  fontsize=12
)

ggsave("Supplemental_fig_bioclim_corr_matrix.png", 
       plot = p, dpi = 300, width = 10, height = 10, units = "in")

###########################################################################
# Find highly correlated variables

cor_matrix <- cor(env_data)

cutoff <- 0.7
vars_to_remove <- findCorrelation(cor_matrix, cutoff = cutoff,names=TRUE)
sort(vars_to_remove)
# BIO3,BIO4, BIO8, BIO15

vars_to_remove

