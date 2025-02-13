# Calculate individual level heterozygosity for inversion regions and genome-wide
# Sophia MacRae Orzechowski
# January 2025



library(ggplot2)
library(dplyr)
setwd("~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/theta")

#########################################################################


theta_data <- read.table("Ben_166_F_bwa_merge_dedup_clip_sort_fixmate_realigned_final.50kb_thetasWindow.gz.pestPG",header=T, sep="\t",comment.char="")
theta_data <- read.table("combined_het_individuals_v2.txt",header=T,sep="\t",comment.char="")
# Convert necessary columns to numeric
theta_data$tH <- as.numeric(theta_data$tH)
theta_data$nSites <- as.numeric(theta_data$nSites)

# Calculate mean heterozygosity per window
#theta_data$Mean_Het <- theta_data$tP / theta_data$nSites
theta_data$Mean_Het <- theta_data$tH / theta_data$nSites

# remove outliers
theta_data <- theta_data[theta_data$Sample_Name!="Binya_366911_M",]

# genome wide heterozygosity
theta_test <- theta_data %>%
  group_by(Sample_Name) %>%
  summarise(
    total_tH = sum(tH, na.rm = TRUE), 
    total_sites = sum(nSites, na.rm = TRUE), 
    overall_het = total_tH / total_sites
  )%>% data.frame()


#########################################################################

## Whole chromosomes
#subset <- theta_data[theta_data$Chr=="scaffold_18_RagTag",]

## outlier regions
#subset <- theta_data[theta_data$Chr=="scaffold_4_RagTag"& theta_data$WinCenter>110000000,]
#subset <- theta_data[theta_data$Chr=="scaffold_6_RagTag"& theta_data$WinCenter>11000000 & theta_data$WinCenter<30000000,]
#subset <- theta_data[theta_data$Chr=="scaffold_8_RagTag"& theta_data$WinCenter>40000000 & theta_data$WinCenter<55000000,]
#subset <- theta_data[theta_data$Chr=="scaffold_8_RagTag"& theta_data$WinCenter<9000000,]
#subset <- theta_data[theta_data$Chr=="scaffold_10_RagTag"& theta_data$WinCenter>13000000,]
#subset <- theta_data[theta_data$Chr=="scaffold_13_RagTag"& theta_data$WinCenter>17000000 & theta_data$WinCenter<22000000,]
#subset <- theta_data[theta_data$Chr=="scaffold_14_RagTag"& theta_data$WinCenter<9000000,]
#subset <- theta_data[theta_data$Chr=="scaffold_15_RagTag"& theta_data$WinCenter<9000000,]
#subset <- theta_data[theta_data$Chr=="scaffold_18_RagTag"& theta_data$WinCenter>4000000& theta_data$WinCenter<9000000,]
#subset <- theta_data[theta_data$Chr=="pri#ptg000045l"& theta_data$WinCenter<10000000,]
#subset <- theta_data[theta_data$Chr=="pri#ptg000015l"& theta_data$WinCenter>10000000,]



## non-outlier regions
#subset <- theta_data[theta_data$Chr=="scaffold_6_RagTag"& theta_data$WinCenter>35000000,]
#subset <- theta_data[theta_data$Chr=="scaffold_8_RagTag"& theta_data$WinCenter<35000000 & theta_data$WinCenter>12000000,]
#subset <- theta_data[theta_data$Chr=="scaffold_13_RagTag"& theta_data$WinCenter<14000000,]
#subset <- theta_data[theta_data$Chr=="scaffold_14_RagTag"& theta_data$WinCenter>11000000,]
#subset <- theta_data[theta_data$Chr=="scaffold_18_RagTag"& theta_data$WinCenter<2000000,]
subset <- theta_data[theta_data$Chr=="pri#ptg000045l"& theta_data$WinCenter>14000000,]
#subset <- theta_data[theta_data$Chr=="pri#ptg000015l"& theta_data$WinCenter<5000000,]


test <- subset %>%
  group_by(Sample_Name) %>%
  summarise(
    total_tH = sum(tH, na.rm = TRUE), 
    total_sites = sum(nSites, na.rm = TRUE), 
    overall_het = total_tH / total_sites
  )%>% data.frame()

colnames(test)[1] <- "Sample"
colnames(theta_test)[1] <- "Sample"
colnames(subset)[1] <- "Sample"

#combine_data <-left_join(test,Chr4.2_combined)
#combine_data <-left_join(test,Chr6_combined)
#combine_data <-left_join(test,Chr8.1_combined)
#combine_data <-left_join(test,Chr8.2_combined)
#combine_data <-left_join(test,Chr10_combined)
#combine_data <-left_join(test,Chr13_combined)
#combine_data <-left_join(theta_test,Chr13_combined)
#combine_data <-left_join(test,Chr14_combined)
#combine_data <-left_join(test,Chr15_combined)
#combine_data <-left_join(test,Chr18_combined)

combine_data <-left_join(test,pri45_combined)
#combine_data <-left_join(test,pri15_combined)

# the rare order for Chr8.1
#combine_data$Cluster_factor <- factor(combine_data$Cluster,levels=c("0","1","2"))

# the usual order for most:
combine_data$Cluster_factor <- factor(combine_data$Cluster,levels=c("0","2","1"))



levels(combine_data$Cluster_factor) <- c("0/0","0/1","1/1")

ggplot(combine_data,aes(Cluster_factor,overall_het))+
  geom_boxplot(notch=FALSE)+
  xlab("")+
  ylab("Heterozygosity")+
 # ylim(0.01,0.014)+ # for pri45
 # ylim(0.0115,0.016)+ # for Chr8
 #ylim(0.0115,0.017)+ # for Chr14
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=16))

#ggplot(combine_data,aes(Cluster_factor,overall_het))+
#  geom_point()



## Create a scatter plot of mean heterozygosity
ggplot(subset, aes(x = WinCenter, y = tH)) +
  geom_point(color = "purple") +
  geom_line(color = "blue", alpha = 0.5) +
  labs(title = "Mean Heterozygosity per Window",
       x = "Window Center",
       y = "Mean Heterozygosity") +
  theme_minimal()+
  theme(legend.position = "none")



combine_subset <-left_join(subset,Chr13_combined)
combine_subset$Cluster <- as.factor(combine_subset$Cluster)
# Create a scatter plot of mean heterozygosity
ggplot(combine_subset, aes(x = WinCenter, y = tH,color=Cluster)) +
  geom_jitter() +
  geom_line() +
  labs(title = "Mean Heterozygosity per Window",
       x = "Window Center",
       y = "Mean Heterozygosity") +
  theme_minimal()+
  theme()+
  scale_color_manual(values=c("blue","red","orange"))
  

