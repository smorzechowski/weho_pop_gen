# Plot theta and Tajima's D per population
# Sophie M Orzechowski
# March 2025

library(ggplot2)
library(reshape2)
library(ggpubr)
library(gghalves)


setwd("~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/theta")


arid <- read.table("arid.50kb_thetasWindow.gz.pestPG",comment.char = "",header=T)
arid$population <- "arid"
arid <- arid[,-1]

temperate <- read.table("temperate.50kb_thetasWindow.gz.pestPG",comment.char = "",header=T)
temperate$population <- "temperate"
temperate <- temperate[,-1]  

theta_data <- rbind(arid,temperate)
  
# Convert necessary columns to numeric
theta_data$tP <- as.numeric(theta_data$tP)
theta_data$nSites <- as.numeric(theta_data$nSites)

# Calculate mean heterozygosity per window
theta_data$Mean_pi <- theta_data$tP / theta_data$nSites


data <- theta_data
data$start <- data$WinCenter -25000

data$Chr <- str_replace(data$Chr,'scaffold','Chr')
data$Chr <- str_remove(data$Chr,'_RagTag')


# genome wide heterozygosity
#theta_test <- theta_data %>%
#  group_by(Sample_Name) %>%
#  summarise(
#    total_tH = sum(tH, na.rm = TRUE), 
#    total_sites = sum(nSites, na.rm = TRUE), 
#    overall_het = total_tH / total_sites
#  )%>% data.frame()

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


data_filt <- data[data$nSites>3000,]


y_top <- max(data_filt$Mean_pi, na.rm = TRUE) * 1.01


mean_pi <- mean(data_filt$Mean_pi)

p <- ggplot(data_filt,aes(start/1000000,Mean_pi,color=population))+
  geom_point(alpha=0.35)+
  geom_segment(data=intervals_df,aes(x=start/1e6,xend=end/1e6),y=y_top,yend=y_top,color='blue',linewidth=2,inherit.aes=FALSE)+
  #ylim(0,0.8)+
  xlab("Chromosome position (Mb)")+
  labs(y = expression(italic(theta[pi])))+
  geom_hline(yintercept = mean_pi, linetype = "dashed", color = "red", size = 1)+
  facet_wrap(.~Chr,scales="free_x")+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="bottom",,
        legend.title=element_blank(),
        legend.text=element_text(size=16))+
  scale_color_manual(values=c("red","black"))



p

ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Pi_all_autosomes_breakpoints.png", 
       plot = p, dpi = 300, width = 14, height = 10, units = "in")

# Just the chromosomes with outliers

data_filt_outliers <- data_filt %>%
  filter(Chr %in% c("Chr_4","Chr_6","Chr_8","Chr_9","Chr_10","Chr_13","Chr_14","Chr_15","Chr_17","Chr_18","pri#ptg000015l","pri#ptg000045l"))

y_top <- max(data_filt_outliers$Mean_pi, na.rm = TRUE) * 1.01

p <- ggplot(data_filt_outliers,aes(start/1000000,Mean_pi,color=population))+
  geom_point(alpha=0.35)+
  geom_segment(data=intervals_df,aes(x=start/1e6,xend=end/1e6),y=y_top,yend=y_top,color='blue',linewidth=2,inherit.aes=FALSE)+
  #ylim(0,0.8)+
  xlab("Chromosome position (Mb)")+
  labs(y = expression(italic(theta[pi])))+
  geom_hline(yintercept = mean_pi, linetype = "dashed", color = "red", size = 1)+
  facet_wrap(.~Chr,scales="free_x")+
  # facet_grid(.~Chr,scales="free_x",space="free_x")+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="bottom",
        legend.title=element_blank(),
        legend.text=element_text(size=16))+
  scale_color_manual(values=c("red","black"))


p

ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Pi_outliers_breakpoints.png", 
       plot = p, dpi = 300, width = 14, height = 10, units = "in")



### Look at Tajima's D in outliers

y_top <- max(data_filt_outliers$Tajima, na.rm = TRUE) * 1.01

p <- ggplot(data_filt_outliers,aes(start/1000000,Tajima,color=population))+
  geom_point(alpha=0.35)+
  geom_segment(data=intervals_df,aes(x=start/1e6,xend=end/1e6),y=y_top,yend=y_top,color='blue',linewidth=2,inherit.aes=FALSE)+
  #ylim(0,0.8)+
  xlab("Chromosome position (Mb)")+
#  stat_smooth(method="loess")+
  labs(y = expression(paste("Tajima's ", italic(D))))+
  geom_hline(yintercept = mean_pi, linetype = "dashed", color = "red", size = 1)+
  facet_wrap(.~Chr,scales="free_x")+
  # facet_grid(.~Chr,scales="free_x",space="free_x")+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="bottom",
        legend.title=element_blank(),
        legend.text=element_text(size=16))+
  scale_color_manual(values=c("red","black"))

p

ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Tajima_outliers_breakpoints.png", 
       plot = p, dpi = 300, width = 14, height = 10, units = "in")


#######################################################################################
### Look at Tajima's D in all chromosomes

mean(data_filt$Tajima[data_filt$in_interval==TRUE & data_filt$population=="arid"])
mean(data_filt$Tajima[data_filt$in_interval==TRUE & data_filt$population=="temperate"])

mean(data_filt$Tajima[data_filt$in_interval==FALSE & data_filt$population=="arid"])
mean(data_filt$Tajima[data_filt$in_interval==FALSE & data_filt$population=="temperate"])

y_top <- max(data_filt$Tajima, na.rm = TRUE) * 1.01

p <- ggplot(data_filt,aes(start/1000000,Tajima,color=population))+
  geom_point(alpha=0.35)+
  geom_segment(data=intervals_df,aes(x=start/1e6,xend=end/1e6),y=y_top,yend=y_top,color='blue',linewidth=2,inherit.aes=FALSE)+
  #ylim(0,0.8)+
  xlab("Chromosome position (Mb)")+
  #ylab("Dxy")+
  labs(y = expression(paste("Tajima's ", italic(D))))+
 # stat_smooth(method="loess")+
 # geom_hline(yintercept = mean_pi, linetype = "dashed", color = "red", size = 1)+
  facet_wrap(.~Chr,scales="free_x")+
  # facet_grid(.~Chr,scales="free_x",space="free_x")+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="bottom",
        legend.title=element_blank(),
        legend.text=element_text(size=16))+
  scale_color_manual(values=c("red","black"))


p

ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Tajima_autosomes_breakpoints.png", 
       plot = p, dpi = 300, width = 14, height = 10, units = "in")



##############################################################################################
# Boxplot of Tajima's D

data_filt$populationf <- as.factor(data_filt$population)

data_filt$Interval[data_filt$in_interval==TRUE]<-"PCRs"
data_filt$Interval[data_filt$in_interval==FALSE]<-"Autosome-wide"

data_filt$Intervalf <- factor(data_filt$Interval)


data_filt <- data_filt[!is.na(data_filt$population),]
data_filt <- data_filt[!is.na(data_filt$in_interval),]


data_filt$populationf <- droplevels(data_filt$populationf)
data_filt$Intervalf   <- droplevels(data_filt$Intervalf)

p <- ggplot(data_filt,aes(Interval,Tajima,fill=population))+
  stat_compare_means(aes(group=population), 
                     method = "t.test",
                     label = "p.signif")+ 
  geom_boxplot(position = position_dodge(width = 0.8),outliers=FALSE)+
  xlab("")+
  labs(y = expression(paste("Tajima's ", italic(D))))+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="bottom",
        legend.title=element_blank(),
        legend.text=element_text(size=16))+
  scale_fill_manual(values=c("red","black"))



# run the t.test manually
t.test(V5 ~ region, data = data[data$region %in% c("Autosome-wide", "neo-Z"), ])
t.test(V5 ~ region, data = data[data$region %in% c("Autosome-wide", "PCRs"), ])
t.test(V5 ~ region, data = data[data$region %in% c("neo-Z", "PCRs"), ])



p

ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Tajima_autosomes_boxplots.png", 
       plot = p, dpi = 300, width = 6, height = 6, units = "in")




table(data_filt$Interval,data_filt$population)

p <- ggplot(data_filt,aes(population,Tajima,fill=Interval))+
  stat_compare_means(aes(group=Interval), 
                     method = "t.test",
                     label = "p.signif")+ 
  geom_boxplot(position = position_dodge(width = 0.8),outliers=FALSE)+
  xlab("")+
  labs(y = expression(italic(Tajima)))+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="bottom")
#scale_color_manual(values=c("black","red"))+




p




####################################################################################
#subset <- data[data$Chr=="pri#ptg000015l",]
#subset <- data[data$Chr=="pri#ptg000045l",]
#subset <- data[data$Chr=="scaffold_4_RagTag",]
subset <- data[data$Chr=="scaffold_8_RagTag",]

subset <- data[data$Chr=="scaffold_15_RagTag",]
subset <- data[data$Chr=="scaffold_18_RagTag",]
subset <- data[data$Chr=="scaffold_10_RagTag",]

# Differences between pops in inversions
subset <- data[data$Chr=="scaffold_6_RagTag",] 
subset <- data[data$Chr=="scaffold_14_RagTag",]

ggplot(subset,aes(WinCenter/1000000,Mean_pi,color=population,fill=population))+
  geom_point()+
  stat_smooth(method="loess")+
  labs(
    x = "Chromosome Position (Mb)",
    y = expression(Theta[pi])  # Adds π as a subscript to Theta
  ) +
    geom_vline(
      data = lines_df[lines_df$Chr=="Chr_8",],
      aes(xintercept = vline_pos / 1e6),
      color = "blue",
      linetype = "dashed",
    )+
  ggtitle("Chromosome 8")+
  theme_minimal()+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=16),
        legend.title=element_text(size=16),
        legend.text=element_text(size=15))




ggplot(theta_data,aes(WinCenter/1000000,Mean_pi,color=population,fill=population))+
  geom_point()+
  facet_wrap(.~Chr,scales="free")+
  stat_smooth(method="loess")+
  xlab("")+
  ylab("Mean Pi")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=16))


#######################################################################################################
# Neo-Z patterns


arid <- read.table("arid_neoZ_males.50kb_thetasWindow.gz.pestPG",comment.char = "",header=T)
arid$population <- "arid"
arid <- arid[,-1]

temperate <- read.table("temperate_neoZ_males.50kb_thetasWindow.gz.pestPG",comment.char = "",header=T)
temperate$population <- "temperate"
temperate <- temperate[,-1]  

theta_data <- rbind(arid,temperate)

# Convert necessary columns to numeric
theta_data$tP <- as.numeric(theta_data$tP)
theta_data$nSites <- as.numeric(theta_data$nSites)

# Calculate mean heterozygosity per window
theta_data$Mean_pi <- theta_data$tP / theta_data$nSites


data <- theta_data
data$start <- data$WinCenter -25000

data$Chr <- str_replace(data$Chr,'scaffold','Chr')
data$Chr <- str_remove(data$Chr,'_RagTag')

data_filt <- data[data$nSites>3000,]


y_top <- max(data_filt$Mean_pi, na.rm = TRUE) * 1.01

mean_pi <- mean(data_filt$Mean_pi)

p <- ggplot(data_filt,aes(start/1000000,Mean_pi,color=population))+
  geom_point(alpha=0.35)+
  xlab("Chromosome position (Mb)")+
  labs(y = expression(italic(theta[pi])))+
  geom_hline(yintercept = mean_pi, linetype = "dashed", color = "red", size = 1)+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="bottom",,
        legend.title=element_blank(),
        legend.text=element_text(size=16))+
  scale_color_manual(values=c("red","black"))



p

ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Pi_neoZ.png", 
       plot = p, dpi = 300, width = 14, height = 10, units = "in")



mean(data_filt$Tajima[data_filt$population=="arid"])
mean(data_filt$Tajima[data_filt$population=="temperate"])


y_top <- max(data_filt$Tajima, na.rm = TRUE) * 1.01

p <- ggplot(data_filt,aes(start/1000000,Tajima,color=population))+
  geom_point(alpha=0.35)+
  xlab("Chromosome position (Mb)")+
  labs(y = expression(paste("Tajima's ", italic(D))))+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="bottom",
        legend.title=element_blank(),
        legend.text=element_text(size=16))+
  scale_color_manual(values=c("red","black"))


p

ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Tajima_neoZ.png", 
       plot = p, dpi = 300, width = 14, height = 10, units = "in")


# Boxplot of Tajima's D


p <- ggplot(data_filt,aes(population,Tajima,fill=population))+
  stat_compare_means(aes(group=population), 
                     method = "t.test",
                     label = "p.signif")+ 
  geom_boxplot(position = position_dodge(width = 0.8),outliers=FALSE)+
  labs(y = expression(paste("Tajima's ", italic(D))))+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="bottom",
        legend.title=element_blank(),
        legend.text=element_text(size=16))+
  scale_fill_manual(values=c("red","black"))


p

ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Tajima_neoZ_boxplot.png", 
       plot = p, dpi = 300, width = 6, height = 6, units = "in")




