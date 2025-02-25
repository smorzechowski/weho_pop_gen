# Plot Fst across the chromosomes
# Sophie C. M. Orzechowski
# February 2025


######################################################################################
library(stringr)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(gghalves)
setwd("~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/fst")
data <- read.table("temperate_arid.pbs.fst.txt",header=F,comment.char = "",skip=1)
data <- read.table("temperate_arid.pbs.fst_50kb.txt",header=F,comment.char = "",skip=1)



######################################################################################

data$Chr <- str_replace(data$V2,'scaffold','Chr')
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
    (data$Chr == chr_i & data$V3 >= start_i & data$V3 <= end_i)
}






data$Chr <- factor(data$Chr,levels=c("Chr_1","Chr_3","Chr_4","Chr_6","Chr_7","Chr_8",
                                     "Chr_9","Chr_10","Chr_11","Chr_12","Chr_13",
                                     "Chr_14","Chr_15","Chr_16","Chr_17","Chr_18",
                                     "Chr_19","Chr_20","Chr_21","Chr_22","Chr_23",
                                     "Chr_24","Chr_25","Chr_26","Chr_27","Chr_28",
                                     "Chr_29","Chr_30","Chr_31","Chr_32","Chr_33",
                                     "Chr_34","Chr_35","Chr_36","Chr_37","Chr_38",
                                     "Chr_39","pri#ptg000015l","pri#ptg000031l",
                                     "pri#ptg000045l","pri#ptg000050l",
                                     "pri#ptg000051l","pri#ptg000085l","pri#ptg000090l"))




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



y_top <- max(data$V5, na.rm = TRUE) * 1.1


######################################################################################

# Inverted chromosomes
data_subset <- data[data$V2=="scaffold_4_RagTag",]
data_subset <- data[data$V2=="scaffold_6_RagTag",]
data_subset <- data[data$V2=="scaffold_8_RagTag",]
data_subset <- data[data$V2=="scaffold_9_RagTag",]
data_subset <- data[data$V2=="scaffold_10_RagTag",]
data_subset <- data[data$V2=="scaffold_13_RagTag",]
data_subset <- data[data$V2=="scaffold_14_RagTag",]
data_subset <- data[data$V2=="scaffold_15_RagTag",]
data_subset <- data[data$V2=="scaffold_17_RagTag",]
data_subset <- data[data$V2=="scaffold_18_RagTag",]
data_subset <- data[data$V2=="pri#ptg000015l",]
data_subset <- data[data$V2=="pri#ptg000045l",]

# Non inverted chromosomes
data_subset <- data[data$V2=="scaffold_1_RagTag",]
data_subset <- data[data$V2=="scaffold_12_RagTag",]
data_subset <- data[data$V2=="scaffold_16_RagTag",]


######################################################################################
# Chromosomes with SVs and contig breakpoints
data_subset <- data[data$V2=="scaffold_8_RagTag",]
intervals_df_subset <- intervals_df[intervals_df$Chr=="Chr_8",]

# Individual plots
p <- ggplot(data_subset,aes(V3/1000000,V5))+
  geom_point()+
  geom_segment(data=intervals_df_subset,aes(x=start/1e6,xend=end/1e6),y=y_top,yend=y_top,color='blue',linewidth=2,inherit.aes=FALSE)+
  ylim(0,1)+
  xlab("Chromosome position (Mb)")+
  #ylab("Fst")+
  labs(y = expression(italic(F)[ST]))+
  geom_vline(xintercept = 1.7, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 10.2, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 15.1, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 58.6, color = "red", linetype = "dashed")+
  theme_minimal()+
  #scale_x_reverse()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="none")
  
p
# Save the plot as a PNG with 300 dpi, for example 10 x 7 inches:
#ggsave("Fst_scaffold_8_RagTag_breakpoints_xflip.svg", plot = p, dpi = 300, width = 10, height = 7, units = "in")
ggsave("Fst_scaffold_8_RagTag_breakpoints.svg", plot = p, dpi = 300, width = 10, height = 7, units = "in")


start_test = 37874992
end_test = 55974910 

data_inv_subset <- data_subset[data_subset$V3>15090000& data_subset$V3<58610000,]
# Individual plots
p <- ggplot(data_inv_subset,aes(V3/1000000,V5))+
  geom_point()+
  geom_segment(x=start_test/1e6,xend=end_test/1e6,y=y_top,yend=y_top,color='blue',linewidth=2)+
  ylim(0,1)+
  xlab("Chromosome position (Mb)")+
  #ylab("Fst")+
  labs(y = expression(italic(F)[ST]))+
 # geom_vline(xintercept = 1.7, color = "red", linetype = "dashed") +
 # geom_vline(xintercept = 10.2, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 15.1, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 58.6, color = "red", linetype = "dashed")+
  theme_minimal()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="none")

# Save the plot as a PNG with 300 dpi, for example 10 x 7 inches:
ggsave("Fst_contig_00047l_breakpoints.svg", plot = p, dpi = 300, width = 10, height = 7, units = "in")


######################################################################################
data_subset <- data[data$V2=="pri#ptg000015l",]
intervals_df_subset <- intervals_df[intervals_df$Chr=="pri#ptg000015l",]

# Individual plots
p <- ggplot(data_subset,aes(V3/1000000,V5))+
  geom_point()+
 geom_segment(data=intervals_df_subset,aes(x=start/1e6,xend=end/1e6),y=y_top,yend=y_top,color='blue',linewidth=2,inherit.aes=FALSE)+
  ylim(0,1)+
  xlab("Chromosome position (Mb)")+
  #ylab("Fst")+
  labs(y = expression(italic(F)[ST]))+
#  geom_vline(xintercept = 1.7, color = "red", linetype = "dashed") +
  theme_minimal()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="none")
p
# Save the plot as a PNG with 300 dpi, for example 10 x 7 inches:
ggsave("Fst_priptg000015l_breakpoints.svg", plot = p, dpi = 300, width = 10, height = 7, units = "in")



######################################################################################
data_subset <- data[data$V2=="scaffold_13_RagTag",]
intervals_df_subset <- intervals_df[intervals_df$Chr=="Chr_13",]

# Individual plots
p <- ggplot(data_subset,aes(V3/1000000,V5))+
  geom_point()+
  geom_segment(data=intervals_df_subset,aes(x=start/1e6,xend=end/1e6),y=y_top,yend=y_top,color='blue',linewidth=2,inherit.aes=FALSE)+
  ylim(0,1)+
  xlab("Chromosome position (Mb)")+
  #ylab("Fst")+
  labs(y = expression(italic(F)[ST]))+
  geom_vline(xintercept = 2.8, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 22.9, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 24.7, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 31.7, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 32.5, color = "red", linetype = "dashed") +
  theme_minimal()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="none")
p
# Save the plot as a PNG with 300 dpi, for example 10 x 7 inches:
ggsave("Fst_scaffold_13_breakpoints.svg", plot = p, dpi = 300, width = 10, height = 7, units = "in")



######################################################################################
data_subset <- data[data$V2=="scaffold_10_RagTag",]
intervals_df_subset <- intervals_df[intervals_df$Chr=="Chr_10",]

# Individual plots
p <- ggplot(data_subset,aes(V3/1000000,V5))+
  geom_point()+
  geom_segment(data=intervals_df_subset,aes(x=start/1e6,xend=end/1e6),y=y_top,yend=y_top,color='blue',linewidth=2,inherit.aes=FALSE)+
  ylim(0,1)+
  xlab("Chromosome position (Mb)")+
  #ylab("Fst")+
  labs(y = expression(italic(F)[ST]))+
  geom_vline(xintercept = 6.3, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 22.3, color = "red", linetype = "dashed") +
  theme_minimal()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="none")
p
# Save the plot as a PNG with 300 dpi, for example 10 x 7 inches:
ggsave("Fst_scaffold_10_breakpoints.svg", plot = p, dpi = 300, width = 10, height = 7, units = "in")


######################################################################################
data_subset <- data[data$V2=="pri#ptg000045l",]
intervals_df_subset <- intervals_df[intervals_df$Chr=="pri#ptg000045l",]

# Individual plots
p <- ggplot(data_subset,aes(V3/1000000,V5))+
  geom_point()+
  geom_segment(data=intervals_df_subset,aes(x=start/1e6,xend=end/1e6),y=y_top,yend=y_top,color='blue',linewidth=2,inherit.aes=FALSE)+
  ylim(0,1)+
  xlab("Chromosome position (Mb)")+
  labs(y = expression(italic(F)[ST]))+
  theme_minimal()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="none")
p
# Save the plot as a PNG with 300 dpi, for example 10 x 7 inches:
ggsave("Fst_priptg000045l_breakpoints.svg", plot = p, dpi = 300, width = 10, height = 7, units = "in")




######################################################################################
data_subset <- data[data$V2=="scaffold_15_RagTag",]
intervals_df_subset <- intervals_df[intervals_df$Chr=="Chr_15",]

# Individual plots
p <- ggplot(data_subset,aes(V3/1000000,V5))+
  geom_point()+
  geom_segment(data=intervals_df_subset,aes(x=start/1e6,xend=end/1e6),y=y_top,yend=y_top,color='blue',linewidth=2,inherit.aes=FALSE)+
  ylim(0,1)+
  xlab("Chromosome position (Mb)")+
  #ylab("Fst")+
  labs(y = expression(italic(F)[ST]))+
  geom_vline(xintercept = 4.5, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 9.2, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 12.3, color = "red", linetype = "dashed") +
  theme_minimal()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="none")
p
# Save the plot as a PNG with 300 dpi, for example 10 x 7 inches:
ggsave("Fst_scaffold_15_breakpoints.svg", plot = p, dpi = 300, width = 10, height = 7, units = "in")


######################################################################################
data_subset <- data[data$V2=="scaffold_6_RagTag",]
intervals_df_subset <- intervals_df[intervals_df$Chr=="Chr_6",]

# Individual plots
p <- ggplot(data_subset,aes(V3/1000000,V5))+
  geom_point()+
  geom_segment(data=intervals_df_subset,aes(x=start/1e6,xend=end/1e6),y=y_top,yend=y_top,color='blue',linewidth=2,inherit.aes=FALSE)+
  ylim(0,1)+
  xlab("Chromosome position (Mb)")+
  #ylab("Fst")+
  labs(y = expression(italic(F)[ST]))+
  geom_vline(xintercept = 8.3, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 16.9, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 31.8, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 72.3, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 74.1, color = "red", linetype = "dashed") +
  theme_minimal()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="none")
p
# Save the plot as a PNG with 300 dpi, for example 10 x 7 inches:
ggsave("Fst_scaffold_6_breakpoints.svg", plot = p, dpi = 300, width = 10, height = 7, units = "in")



######################################################################################
# All chromosomes

ggplot(data,aes(V3/1000000,V5,color=in_interval))+
  geom_point()+
  geom_segment(data=intervals_df,aes(x=start/1e6,xend=end/1e6),y=y_top,yend=y_top,color='blue',linewidth=2,inherit.aes=FALSE)+
  ylim(0,0.8)+
  xlab("Chromosome position (Mb)")+
  #ylab("Fst")+
  labs(y = expression(italic(F)[ST]))+
  facet_wrap(.~Chr,scales="free_x")+
 # facet_grid(.~Chr,scales="free_x",space="free_x")+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="none")+
  scale_color_manual(values=c("grey","darkblue"))

######################################################################################
# box plot of Fst inside and outside inversions
ggplot(data,aes(in_interval,V5))+
  geom_boxplot()+
  ylim(0,1)+
  xlab("")+
  labs(y = expression(italic(F)[ST]))+
  stat_compare_means(comparisons = list(c("TRUE", "FALSE")), 
                     method = "t.test") +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        legend.position="none")+
  scale_x_discrete(labels = c("TRUE" = "Inversions", "FALSE" = "Genome-wide"))

######################################################################################
# Histogram of Fst inside and outside inversions
ggplot(data,aes(V5,fill=in_interval,color=in_interval))+
  geom_histogram(bins=100,alpha=.3,position='identity')+
  ylab("Count (50 kb windows)")+
  labs(x = expression(italic(F)[ST]))+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        legend.position="none")





