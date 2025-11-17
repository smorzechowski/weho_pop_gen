#Load the RcppCNPy package
library(RcppCNPy)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(patchwork)
library(pals)

setwd("~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/pcangsd")
samples<-read.table('~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/ngsRelate/allbamspath.txt',sep="/")
samples$samples <- str_remove_all(samples$V12,'_bwa_merge_dedup_clip_sort_fixmate_realigned_final.bam|_bwa_dedup_clip_sort_fixmate_realigned_final.bam')

samples_excl_1st <- samples
samples_excl_1st <- samples[samples$samples!="Berth_162_F" &
                            samples$samples!="Mallee_313_U" &
                            samples$samples!="Moon_365024_M" &
                            samples$samples!="Nom_100_M" &
                            samples$samples!="Pill_110_M" &
                            samples$samples!="Pill_111_M" &
                            samples$samples!="Walch_154_M",]

#### Autosomal structure, excluding inversions and first order relatives and only 39 scaffolds included
#C <- as.matrix(read.table("Nleu_autos_rm_inversions_filtpcangsd.cov")) # Reads estimated covariance matrix
#C <- as.matrix(read.table("Nleu_autos_rm_inversions_filt_excl_1st_thin_pcangsd_thinned.cov")) # Reads estimated covariance matrix
C <- as.matrix(read.table("Nleu_autos_rm_inversions_filt_excl_1st_thin_39scaf_pcangsd_thinned.cov")) # Reads estimated covariance matrix

Apop_excl_1st <- c("Weddin","Weddin","Weddin","Weddin","Weddin","Weddin","Weddin","Binya","Binya","Binya","Binya","Binya","Binya",
                   "Bog","Bog","Weddin","Gund","Gund","Gund","Ing","Ing","Ing","Mallee","Mallee","Mallee","Mallee","Mallee","Mallee","Mallee",
                   "Moon","Moon","Moon","Moon","Mull","Mull","Mull","Mull","Mull","Mull","Nom","Nom","Nom","Nom","Nom",
                   "Nom","Nom","Pill","Pill","Pill","Pill","Pill","Pill","Reedy","Reedy","Reedy","Reedy","Talla",
                   "Walch","Walch","Walch","Walch","Walch","Walch","Walch","Zost","Zost","Zost")

#Apop <- c("Weddin","Weddin","Weddin","Weddin","Weddin","Weddin","Weddin","Weddin","Binya","Binya","Binya","Binya","Binya","Binya",
#                   "Bog","Bog","Weddin","Gund","Gund","Gund","Ing","Ing","Ing","Mall","Mall","Mall","Mall","Mall","Mall","Mall","Mall",
#                   "Moon","Moon","Moon","Moon","Moon","Mull","Mull","Mull","Mull","Mull","Mull","Nom","Nom","Nom","Nom","Nom","Nom",
#                   "Nom","Nom","Pill","Pill","Pill","Pill","Pill","Pill","Pill","Pill","Reedy","Reedy","Reedy","Reedy","Talla",
#                   "Walch","Walch","Walch","Walch","Walch","Walch","Walch","Walch","Zost","Zost","Zost")


mme.pca <- eigen(C) #perform the pca using the eigen function. 

eigenvectors = mme.pca$vectors #extract eigenvectors 
pca.vectors = as_tibble(cbind(Apop_excl_1st, data.frame(eigenvectors))) #combine with our population assignments

# subset to get 65 individuals in RDA analysis
pca.vectors.samples = as_tibble(cbind(samples_excl_1st$samples,pca.vectors)) #combine with our sample assignments
colnames(pca.vectors.samples)[1] <- "Sample"
remove <- c("Binya_366911_M","Walch_120_F")
pca.vectors.samples.subset <- pca.vectors.samples[!is.element(pca.vectors.samples$Sample,remove),]
PC1 <- pca.vectors.samples.subset$X1
write.table(PC1,"Nleu_PC1_65ind_autos_excl_inv.txt",row.names = FALSE,col.names = FALSE,quote = FALSE)


pca.eigenval.sum = sum(mme.pca$values) #sum of eigenvalues
varPC1 <- (mme.pca$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1
varPC2 <- (mme.pca$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2
varPC3 <- (mme.pca$values[3]/pca.eigenval.sum)*100 #Variance explained by PC3
varPC4 <- (mme.pca$values[4]/pca.eigenval.sum)*100 #Variance explained by PC4
varPC5 <- (mme.pca$values[5]/pca.eigenval.sum)*100 #Variance explained by PC4
varPC6 <- (mme.pca$values[6]/pca.eigenval.sum)*100 #Variance explained by PC4
varPC67 <- (mme.pca$values[67]/pca.eigenval.sum)*100 #Variance explained by PC4

length(mme.pca$values)

varPC1
varPC2
varPC3
varPC4
varPC5
varPC6

varPC67

Apop_excl_1st_f <- factor(Apop_excl_1st)

colors <- unname(alphabet()[1:14])
colors <- unname(glasbey()[1:14])

pca_exc <- ggplot(data = pca.vectors, aes(x = X1, y = X2, colour = Apop_excl_1st_f)) +
  geom_point(size = 3) +
 ylim(-0.4,0.4)+
 xlim(-0.3,0.3)+
 # geom_text(aes(label = samples_excl_1st$samples), vjust = -1, size = 3,show.legend=FALSE) +
  #geom_text_repel(aes(label = pop), size = 3,max.overlaps = 10) +
  labs(title = "A.",# Autosomal population structure excluding inversions",
       x = paste0("PC1 (", round(varPC1, 1), "%)"),
       y = paste0("PC2 (", round(varPC2, 1), "%)")) +
  theme_minimal()+
  scale_x_reverse(limits = c(0.3,-0.3), breaks = c(-0.2,-0.1,0,0.1,0.2))+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=16),
        plot.title=element_text(face="bold",size=16),
        legend.title=element_blank(),
        legend.text=element_text(size=15),
        legend.position="none")+
  scale_color_manual(values=colors)



print(pca_exc)



#### Autosomal structure, including inversions, excluding 1st order, including only 39 scaffolds
#C <- as.matrix(read.table("Nleu_autos_ANGSD_PCApcangsd.cov")) # Reads estimated covariance matrix
#C <- as.matrix(read.table("Nleu_autos_rm_filt_excl_1st_thin_pcangsd_thinned.cov")) # Reads estimated covariance matrix
C <- as.matrix(read.table("Nleu_autos_rm_filt_excl_1st_thin_39scaf_pcangsd_thinned.cov")) # Reads estimated covariance matrix

#Apop <- c("Weddin","Weddin","Weddin","Weddin","Weddin","Weddin","Weddin","Weddin","Binya","Binya","Binya","Binya","Binya","Binya",
#          "Bog","Bog","Weddin","Gund","Gund","Gund","Ing","Ing","Ing","Mallee","Mallee","Mallee","Mallee","Mallee","Mallee","Mallee","Mallee",
#          "Moon","Moon","Moon","Moon","Moon","Mull","Mull","Mull","Mull","Mull","Mull","Nom","Nom","Nom","Nom","Nom","Nom",
#          "Nom","Nom","Pill","Pill","Pill","Pill","Pill","Pill","Pill","Pill","Reedy","Reedy","Reedy","Reedy","Talla",
#          "Walch","Walch","Walch","Walch","Walch","Walch","Walch","Walch","Zost","Zost","Zost")


Apop_excl_1st <- c("Weddin","Weddin","Weddin","Weddin","Weddin","Weddin","Weddin","Binya","Binya","Binya","Binya","Binya","Binya",
                   "Bog","Bog","Weddin","Gund","Gund","Gund","Ing","Ing","Ing","Mallee","Mallee","Mallee","Mallee","Mallee","Mallee","Mallee",
                   "Moon","Moon","Moon","Moon","Mull","Mull","Mull","Mull","Mull","Mull","Nom","Nom","Nom","Nom","Nom",
                   "Nom","Nom","Pill","Pill","Pill","Pill","Pill","Pill","Reedy","Reedy","Reedy","Reedy","Talla",
                   "Walch","Walch","Walch","Walch","Walch","Walch","Walch","Zost","Zost","Zost")

mme.pca <- eigen(C) #perform the pca using the eigen function. 

eigenvectors = mme.pca$vectors #extract eigenvectors 
pca.vectors = as_tibble(cbind(Apop_excl_1st, data.frame(eigenvectors))) #combine with our population assignments


pca.eigenval.sum = sum(mme.pca$values) #sum of eigenvalues
varPC1 <- (mme.pca$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1
varPC2 <- (mme.pca$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2
varPC3 <- (mme.pca$values[3]/pca.eigenval.sum)*100 #Variance explained by PC3
varPC4 <- (mme.pca$values[4]/pca.eigenval.sum)*100 #Variance explained by PC4
varPC5 <- (mme.pca$values[5]/pca.eigenval.sum)*100 #Variance explained by PC4
varPC6 <- (mme.pca$values[6]/pca.eigenval.sum)*100 #Variance explained by PC4

varPC1
varPC2
varPC3
varPC4
varPC5
varPC6


Apop_excl_1st_f <- factor(Apop_excl_1st)

#colors <- unname(alphabet()[1:14])
colors <- unname(glasbey()[1:14])


pca_inc <- ggplot(data = pca.vectors, aes(x = X1, y = X2, colour = Apop_excl_1st)) +
  geom_point(size = 3) +
  ylim(-0.4,0.4)+
  xlim(-0.3,0.3)+
#  geom_text(aes(label = samples_excl_1st$samples), vjust = -1, size = 3,show.legend=FALSE) +
  #geom_text_repel(aes(label = pop), size = 3,max.overlaps = 10) +
  labs(title = "B.",# Autosomal population structure including inversions",
       x = paste0("PC1 (", round(varPC1, 1), "%)"),
       y = paste0("PC2 (", round(varPC2, 1), "%)")) +
  theme_minimal()+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=16),
        plot.title=element_text(face="bold",size=16),
        legend.title=element_blank(),
        legend.text=element_text(size=15),
        legend.position="none")+
  scale_color_manual(values=colors)+
  scale_x_continuous(limits = c(-0.3,0.3), breaks = c(-0.2,-0.1,0,0.1,0.2))

print(pca_inc)


#p <- pca_inc/pca_exc
p <- pca_exc/pca_inc

p
ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Figure3_Pcangsd_autos_PCA_with_without_inversions_flipped.png", 
       plot = p, dpi = 300, width = 6, height = 8, units = "in")

ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Figure3_Pcangsd_autos_PCA_with_without_inversions.png", 
       plot = p, dpi = 300, width = 6, height = 8, units = "in")





#####################################################################################################

#### New PAR structure, excluding first order relatives
#C <- as.matrix(read.table("Nleu_newPAR_pcangsd_unpruned.cov")) # Reads estimated covariance matrix
#C <- as.matrix(read.table("Nleu_newPAR_rm_filt_excl_1st_order_pcangsd_unpruned.cov")) # Reads estimated covariance matrix
C <- as.matrix(read.table("Nleu_newPAR_rm_filt_excl_1st_order_thin_updated_pcangsd_thinned.cov")) # Reads estimated covariance matrix

#Apop <- c("Weddin","Weddin","Weddin","Weddin","Weddin","Weddin","Weddin","Weddin","Binya","Binya","Binya","Binya","Binya","Binya",
#          "Bog","Bog","Weddin","Gund","Gund","Gund","Ing","Ing","Ing","Mallee","Mallee","Mallee","Mallee","Mallee","Mallee","Mallee","Mallee",
#          "Moon","Moon","Moon","Moon","Moon","Mull","Mull","Mull","Mull","Mull","Mull","Nom","Nom","Nom","Nom","Nom","Nom",
#          "Nom","Nom","Pill","Pill","Pill","Pill","Pill","Pill","Pill","Pill","Reedy","Reedy","Reedy","Reedy","Talla",
#          "Walch","Walch","Walch","Walch","Walch","Walch","Walch","Walch","Zost","Zost","Zost")

Apop_excl_1st <- c("Weddin","Weddin","Weddin","Weddin","Weddin","Weddin","Weddin","Binya","Binya","Binya","Binya","Binya","Binya",
                   "Bog","Bog","Weddin","Gund","Gund","Gund","Ing","Ing","Ing","Mallee","Mallee","Mallee","Mallee","Mallee","Mallee","Mallee",
                   "Moon","Moon","Moon","Moon","Mull","Mull","Mull","Mull","Mull","Mull","Nom","Nom","Nom","Nom","Nom",
                   "Nom","Nom","Pill","Pill","Pill","Pill","Pill","Pill","Reedy","Reedy","Reedy","Reedy","Talla",
                   "Walch","Walch","Walch","Walch","Walch","Walch","Walch","Zost","Zost","Zost")

mme.pca <- eigen(C) #perform the pca using the eigen function. 

eigenvectors = mme.pca$vectors #extract eigenvectors 
pca.vectors = as_tibble(cbind(Apop_excl_1st, data.frame(eigenvectors))) #combine with our population assignments


pca.eigenval.sum = sum(mme.pca$values) #sum of eigenvalues
varPC1 <- (mme.pca$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1
varPC2 <- (mme.pca$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2
varPC3 <- (mme.pca$values[3]/pca.eigenval.sum)*100 #Variance explained by PC3
varPC4 <- (mme.pca$values[4]/pca.eigenval.sum)*100 #Variance explained by PC4

varPC1
varPC2
varPC3
varPC4


Apop_excl_1st_f <- factor(Apop_excl_1st)

colors <- unname(alphabet()[1:14])
colors <- unname(glasbey()[1:14])


pca_par <- ggplot(data = pca.vectors, aes(x = X1, y = X2, colour = Apop_excl_1st_f)) +
  geom_point(size = 3) +
  ylim(-0.7,0.7)+
  xlim(-0.5,0.5)+
  #geom_text(aes(label = samples_excl_1st$samples), vjust = -1, size = 3,show.legend=FALSE) +
  #geom_text_repel(aes(label = pop), size = 3,max.overlaps = 10) +
  labs(title ="C.", "New PAR population structure",
       x = paste0("PC1 (", round(varPC1, 1), "%)"),
       y = paste0("PC2 (", round(varPC2, 1), "%)")) +
  theme_minimal()+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=16),
        plot.title=element_text(face="bold",size=16),
        legend.title=element_blank(),
        legend.text=element_text(size=15),
        legend.position="none")+
  scale_color_manual(values=colors)#+
 # scale_x_continuous(limits = c(-0.3,0.3), breaks = c(-0.2,-0.1,0,0.1,0.2))

print(pca_par)

########################################################################################################

#### neo-Z structure in males, excluding 1st order relatives and thinned
#C <- as.matrix(read.table("Nleu_neoZ_males_pcangsd_unpruned.cov")) # Reads estimated covariance matrix
C <- as.matrix(read.table("Nleu_neoZ_males_rm_filt_excl_1st_order_thin_pcangsd_thinned.cov")) # Reads estimated covariance matrix

Mpop_excl_1st <- c("Nom","Weddin","Weddin","Weddin","Binya","Binya","Binya","Binya","Bog","Weddin","Gund","Gund","Gund","Mallee","Mallee",
          "Mallee","Mallee","Mallee","Moon","Moon","Moon","Mull","Mull","Mull","Mull","Mull","Nom","Nom","Nom","Nom",
          "Pill","Pill","Pill","Pill","Pill","Reedy","Reedy","Reedy","Walch","Walch","Walch","Walch","Zost","Zost")


Msam_excl_1st <- c("Nom-93","Weddin","Weddin","Weddin","Binya","Binya","Binya","Binya","Bog","Weddin","Gund","Gund","Gund","Mallee","Mallee",
          "Mallee","Mallee","Mallee","Moon","Moon","Moon","Mull-107","Mull-918","Mull-919","Mull-920","Mull-921","Nom-96","Nom-100","Nom-103","Nom-105",
          "Pill","Pill","Pill","Pill","Pill","Reedy","Reedy","Reedy","Walch","Walch","Walch","Walch","Zost","Zost")

mme.pca <- eigen(C) #perform the pca using the eigen function. 

eigenvectors = mme.pca$vectors #extract eigenvectors 
pca.vectors = as_tibble(cbind(Mpop_excl_1st, data.frame(eigenvectors))) #combine with our population assignments

pca.eigenval.sum = sum(mme.pca$values) #sum of eigenvalues
varPC1 <- (mme.pca$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1
varPC2 <- (mme.pca$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2
varPC3 <- (mme.pca$values[3]/pca.eigenval.sum)*100 #Variance explained by PC3
varPC4 <- (mme.pca$values[4]/pca.eigenval.sum)*100 #Variance explained by PC4

varPC1
varPC2
varPC3
varPC4

Mpop_excl_1st_f <- factor(Mpop_excl_1st)

colors <- unname(alphabet()[1:14])
colors <- unname(glasbey()[1:14])

Mcolors <- c("#0000FF","#FF0000","#00FF00","#FF00B6","#005300","#FFD300","#009FFF","#9A4D42","#00FFBE","#1F9698","#FFACFD","#B1CC71")


pca_neoZ <- ggplot(data = pca.vectors, aes(x = X1, y = X2, colour = Mpop_excl_1st_f)) +
  geom_point(size = 3) +
  ylim(-0.7,0.7)+
 # xlim(-0.5,0.5)+
 # geom_text(aes(label = Msam_excl_1st), vjust = -1, size = 3,show.legend=FALSE) +
  #geom_text_repel(aes(label = pop), size = 3,max.overlaps = 10) +
  labs(title = "D.", #neo-Z population structure in males",
       x = paste0("PC1 (", round(varPC1, 1), "%)"),
       y = paste0("PC2 (", round(varPC2, 1), "%)")) +
  theme_minimal()+
  theme(axis.text=element_text(size=15),
        plot.title=element_text(face="bold",size=16),
        axis.title=element_text(size=16),
        legend.title=element_blank(),
        legend.text=element_text(size=15),
        legend.position="none")+
  scale_color_manual(values=Mcolors)+
  scale_x_reverse(limits = c(0.5,-0.5), breaks = c(-0.5,-0.25,0,0.25,0.5))

print(pca_neoZ)

p <- (pca_par/pca_neoZ)

p
ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Figure4_Pcangsd_neoZ_newPAR.png", 
       plot = p, dpi = 300, width = 6, height = 8, units = "in")


# Final figure: all combined into one - can't combine legends with plot_layout(guides = "collect") because they are different

#final <- (pca_exc + pca_inc) / (pca_par + pca_neoZ) +
#  plot_layout(guides = "collect")+
#  theme(
#    legend.position = "top")
#    legend.justification = "center",
#    plot.margin = margin(5, 40, 5, 5)   # add space on the right
#  )

final <- (pca_exc + pca_inc) / (pca_par + pca_neoZ)+
plot_layout(guides = "collect") +
  theme(legend.position = "bottom")

final

#######################
# Create dummy data for legend
sites <- c("Binya","Bog","Gund","Ing","Mallee","Moon","Mull","Nom","Pill","Reedy","Talla","Walch","Weddin","Zost")
cols <- unname(glasbey()[1:14])        # 12 colors from a palette
#cols <- c(cols, "#999999", "#000000")                 # add 2 more manually

legend_df <- data.frame(site = factor(sites, levels = sites))

legend_df$site
legend_df$site <- factor(legend_df$site, levels = c("Mallee","Binya","Ing","Nom","Weddin","Zost","Gund","Pill","Mull","Reedy","Talla","Moon","Bog","Walch"))
cols_reordered <- cols[match(levels(legend_df$site), sites)]
# Make an empty plot that only draws points (for legend)
dummy_legend <- ggplot(legend_df, aes(x = 1, y = 1, color = site)) +
  geom_point(size = 4) +
  scale_color_manual(values = cols_reordered, name = "Location") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold",size=15),
    legend.key.size = unit(0.6, "lines"),
    legend.text = element_text(size = 14)
  )

dummy_legend


components <- ggplotGrob(dummy_legend)
custom_legend <- gtable::gtable_filter(components, "guide-box")

final_legend <- plot_grid(final,custom_legend,ncol = 2,rel_widths = c(1,0.15))
final_legend
# WEHE pop gen folder
ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Figure2_Pcangsd_autos_neoZ_newPAR.png", 
       plot = final, dpi = 300, width = 10, height = 8, units = "in")

# Manuscript folder
ggsave("C:/Users/sophi/Documents/PhD research/Manuscripts/WEHE pop gen manuscript/figures/Figure2_Pcangsd_autos_neoZ_newPAR.png", 
       plot = final, dpi = 300, width = 10, height = 8, units = "in")


## Combine PCA with admixture plots

## Create plot_admix and plot_admix_legend with PCangsd_admix.R

combined_final <- (final / plot_admix_legend)
combined_final

combined_final <- (final_legend / plot_admix)
combined_final

# Manuscript folder
ggsave("C:/Users/sophi/Documents/PhD research/Manuscripts/WEHE pop gen manuscript/figures/Figure2_Pcangsd_autos_neoZ_newPAR.png", 
       plot = combined_final, dpi = 300, width = 8.5, height = 11, units = "in")
