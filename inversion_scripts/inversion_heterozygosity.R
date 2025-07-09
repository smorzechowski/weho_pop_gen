# Calculate individual level heterozygosity for inversion regions and genome-wide
# Sophia MacRae Orzechowski
# January 2025



library(ggplot2)
library(dplyr)
library(stringr)

source("~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE R scripts/inversion_clusters_add_latlong.R")
setwd("~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/theta")

#########################################################################


theta_data <- read.table("Ben_166_F_bwa_merge_dedup_clip_sort_fixmate_realigned_final.50kb_thetasWindow.gz.pestPG",header=T, sep="\t",comment.char="")
theta_data <- read.table("combined_het_individuals_v2.txt",header=T,sep="\t",comment.char="")



# PCR regions with alternate patterns

#Chr9_combined
#Chr8.1_combined

# Additional PCR on Chr4
#Chr4.1_combined

# Change to Chr_8 so that the join below works better
Chr4.2_combined$Chr <- "Chr_4"


# Manually fix levels for Chr9
Chr9_combined$Cluster_fix <- Chr9_combined$Cluster
Chr9_combined$Cluster_fix[Chr9_combined$Cluster==0] <- 2
Chr9_combined$Cluster_fix[Chr9_combined$Cluster==2] <- 0
Chr9_combined$Cluster <- NULL
colnames(Chr9_combined)[6] <- "Cluster"

# Manually fix levels for Chr8.1
Chr8.1_combined$Cluster_fix <- Chr8.1_combined$Cluster
Chr8.1_combined$Cluster_fix[Chr8.1_combined$Cluster==1] <- 2
Chr8.1_combined$Cluster_fix[Chr8.1_combined$Cluster==2] <- 1
Chr8.1_combined$Cluster <- NULL
colnames(Chr8.1_combined)[6] <- "Cluster"


cluster_data <- rbind(Chr4.2_combined,Chr6_combined,Chr8.1_combined,Chr8.2_combined, Chr9_combined,Chr10_combined,
                      Chr13_combined,Chr14_combined,Chr15_combined,Chr17_combined,Chr18_combined,pri15_combined,pri45_combined)

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


data <- theta_data
data$start <- data$WinCenter -25000

data$Chr <- str_replace(data$Chr,'scaffold','Chr')
data$Chr <- str_remove(data$Chr,'_RagTag')


# Create data frame for vertical lines 
lines_df <- data.frame(
  Chr = c(
    # "Chr_4", "Chr_4",
    "Chr_4", "Chr_4",
    "Chr_6", "Chr_6",
    # "Chr_6", "Chr_6",
    "Chr_8", "Chr_8",
    "Chr_8", "Chr_8",
    "Chr_9", "Chr_9",
    # "Chr_9", "Chr_9",
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
    55974910, 37874992,
    # Chr_9_:16374481-25249745
    1631342, 18339159,
    # Chr_9_:16374481-18339159
    #16374481, 18339159,
    # pri#ptg000015l:6864979-17909547
    6864979, 17909547,
    # Chr_10_:9675047-22345924
    9675047, 22345924,
    # Chr_13_:15724920-25076809
    # Chr_13_:14500000-25076809
    # Updated to 23500000
    15000000, 23500000,
    # Chr_14_:1-101749905
    1, 10174990,
    # Chr_15_:1-877497
    1, 8774970,
    # Chr_17_:12914902-19617658
    12708285, 19617658,
    # Chr_18_:1865069-10487981
    # updated this!!
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



data$in_interval <- FALSE

data$Chr <- as.character(data$Chr)
intervals_df$Chr <- as.character(intervals_df$Chr)

# For each interval, mark those points as in_interval
for(i in seq_len(nrow(intervals_df))) {
  chr_i   <- intervals_df$Chr[i]
  start_i <- intervals_df$start[i]
  end_i   <- intervals_df$end[i]
  
  data$in_interval <- data$in_interval | 
    (data$Chr == chr_i & data$start >= start_i & data$start <= end_i)
}


#table(data$in_interval)

mean_het <- mean(data$Mean_pi)

data$Chr[data$Chr=="Chr_8" & data$start<10724985] <- "Chr_8.2"
data$Chr[data$Chr=="Chr_8" & data$start>37874992 & data$start<55974910] <- "Chr_8.1"


test <- data %>%
  group_by(Sample_Name,Chr,in_interval) %>%
  summarise(
    total_tH = sum(tH, na.rm = TRUE), 
    total_sites = sum(nSites, na.rm = TRUE), 
    overall_het = total_tH / total_sites
  )%>% data.frame()


colnames(test)[1] <- "Sample"

combine_data <-left_join(test,cluster_data)


# the rare order for Chr8.1
#combine_data$Cluster_factor <- factor(combine_data$Cluster,levels=c("0","1","2"))

# the usual order for most:
combine_data$Cluster_factor <- factor(combine_data$Cluster,levels=c("0","2","1"))

levels(combine_data$Cluster_factor) <- c("0/0","0/1","1/1")


data_filt <- combine_data[combine_data$in_interval==TRUE,]
#data_filt <- data_filt[data_filt$Chr!="Chr_4"&data_filt$Chr!="Chr_8",]
#data_filt <- data_filt[data_filt$Chr!="Chr_9",]
#data_filt <- data_filt[data_filt$Chr!="Chr_8",]

data_filt$Chr <- factor(data_filt$Chr,levels=c("Chr_1","Chr_3","Chr_4","Chr_6","Chr_7","Chr_8.1","Chr_8.2",
                                     "Chr_9","Chr_10","Chr_11","Chr_12","Chr_13",
                                     "Chr_14","Chr_15","Chr_16","Chr_17","Chr_18",
                                     "Chr_19","Chr_20","Chr_21","Chr_22","Chr_23",
                                     "Chr_24","Chr_25","Chr_26","Chr_27","Chr_28",
                                     "Chr_29","Chr_30","Chr_31","Chr_32","Chr_33",
                                     "Chr_34","Chr_35","Chr_36","Chr_37","Chr_38",
                                     "Chr_39","pri#ptg000015l","pri#ptg000031l",
                                     "pri#ptg000045l","pri#ptg000050l",
                                     "pri#ptg000051l","pri#ptg000085l","pri#ptg000090l"))


p <- ggplot(data_filt,aes(Cluster_factor,overall_het))+
  geom_boxplot()+
  xlab("")+
  labs(y = expression(italic(H[o])))+
  facet_wrap(.~Chr,scales="free_x")+
  theme_minimal()+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="bottom",,
        legend.title=element_blank(),
        legend.text=element_text(size=16))+
  scale_color_manual(values=c("red","black"))

p


ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/heterozygosity_clusters_incl_chr8.png", 
       plot = p, dpi = 350, width = 10, height = 10, units = "in")





#########################################################################

## Whole chromosomes
#subset <- theta_data[theta_data$Chr=="scaffold_18_RagTag",]

## outlier regions
# NOTE these are very ROUGH subsets compared to the in_interval variable I made above -- this is why plots might not match exactly
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
  

