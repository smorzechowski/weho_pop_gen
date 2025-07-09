# Plot theta and Tajima's D per population
# Sophie M Orzechowski
# March 2025

library(ggplot2)
library(reshape2)
library(ggpubr)
library(gghalves)
library(stringr)
library(dplyr)


setwd("~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/theta")


arid <- read.table("arid.50kb_thetasWindow.gz.pestPG",comment.char = "",header=T)
arid$population <- "arid"
arid <- arid[,-1]

temperate <- read.table("temperate.50kb_thetasWindow.gz.pestPG",comment.char = "",header=T)
temperate$population <- "temperate"
temperate <- temperate[,-1]  


aridZ <- read.table("arid_neoZ_males.50kb_thetasWindow.gz.pestPG",comment.char = "",header=T)
aridZ$population <- "arid"
aridZ <- aridZ[,-1]

temperateZ <- read.table("temperate_neoZ_males.50kb_thetasWindow.gz.pestPG",comment.char = "",header=T)
temperateZ$population <- "temperate"
temperateZ <- temperateZ[,-1]  

theta_data <- rbind(arid,aridZ,temperate,temperateZ)
theta_data$start <- theta_data$WinCenter -25000


# Convert necessary columns to numeric
theta_data$tP <- as.numeric(theta_data$tP)
theta_data$nSites <- as.numeric(theta_data$nSites)

# Calculate mean heterozygosity per window
theta_data$Mean_pi <- theta_data$tP / theta_data$nSites


data <- theta_data
data$start <- data$WinCenter -25000
data$Chr[data$Chr=="scaffold_2_RagTag"]<-"Neo_Z"
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



data$Chr <- factor(data$Chr,levels=c("Neo_Z","Chr_1","Chr_3","Chr_4","Chr_6","Chr_7","Chr_8",
                                     "Chr_9","Chr_10","Chr_11","Chr_12","Chr_13",
                                     "Chr_14","Chr_15","Chr_16","Chr_17","Chr_18",
                                     "Chr_19","Chr_20","Chr_21","Chr_22","Chr_23",
                                     "Chr_24","Chr_25","Chr_26","Chr_27","Chr_28",
                                     "Chr_29","Chr_30","Chr_31","Chr_32","Chr_33",
                                     "Chr_34","Chr_35","Chr_36","Chr_37","Chr_38",
                                     "Chr_39","pri#ptg000015l","pri#ptg000031l",
                                     "pri#ptg000045l","pri#ptg000050l",
                                     "pri#ptg000051l","pri#ptg000085l","pri#ptg000090l"))




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


data_filt <- data[data$nSites>3000,]





mean_pi <- mean(data_filt$Mean_pi)

#########################################################################
# Create dummy legend

# Create the dummy legend data frame
df_legend <- data.frame(
  x = c(1, 1, 1),      
  y = c(1, 2, 3),      
  segment_type = factor(c("MORs", "arid","temperate"), levels = c("MORs", "arid","temperate"))
)

# Generate the dummy plot with both segments & the SNP circle
dummy_legend_plot <- ggplot(df_legend, aes(x, y, fill = segment_type, color = segment_type, shape = segment_type)) +
  # Rectangular legend item for "PCRs"
  geom_tile(data = df_legend[df_legend$segment_type == "MORs", ], width = 0.7, height = 0.3) +
  
  # Large red circle legend item for "dummy SNP"
  geom_point(data = df_legend[df_legend$segment_type != "MORs", ], 
             size = 1,alpha=0.35) +
  
  # Define manual scales for color, fill, and shape
  scale_fill_manual(
    name = "Legend",  # Legend title
    values = c("MORs" = "blue", "arid" = "red", "temperate" = "black")  # Both must be defined here
  ) +
  scale_color_manual(
    name = "Legend",
    values = c("MORs" = "blue", "arid" = "red", "temperate" = "black")  # Colors must be consistent
  ) +
  scale_shape_manual(
    name = "Legend",
    values = c("MORs" = NA, "arid" = 19, "temperate" = 19)  # Shape 21 is a circle; 16 is a filled circle
  ) +
  
  # Force ggplot to combine them into a single legend
  guides(
    fill = guide_legend(override.aes = list(shape = c(NA, 19, 19), color = c("blue", "red","black"), size = c(NA, 1.5, 1.5),alpha=c(NA,0.35,0.35))),
    color = "none",  # Hide redundant color legend
    shape = "none"   # Hide redundant shape legend
  ) +
  
  theme_void() +
  theme(legend.position = "top", legend.title = element_blank())

# Display the legend
print(dummy_legend_plot)

# Extract the legend
components <- ggplotGrob(dummy_legend_plot)

# Print all grobs in the plot to see where the legend is
print(components)
custom_legend <- gtable::gtable_filter(components, "guide-box")

# Print the extracted legend (This is now a gtable, NOT a full ggplot)
grid::grid.draw(custom_legend)

y_top <- max(data_filt$Mean_pi, na.rm = TRUE) * 1.01

p <- ggplot(data_filt,aes(start/1000000,Mean_pi,color=population))+
  geom_point(alpha=0.35)+
  geom_segment(data=intervals_df,aes(x=start/1e6,xend=end/1e6),y=y_top,yend=y_top,color='blue',linewidth=2,inherit.aes=FALSE)+
  xlab("Chromosome position (Mb)")+
  labs(y = expression(italic(theta[pi])))+
  geom_hline(yintercept = mean_pi, linetype = "dashed", color = "red", size = 1)+
  facet_wrap(.~Chr,scales="free_x")+
  theme_minimal()+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="none",,
        legend.title=element_blank(),
        legend.text=element_text(size=16),
        panel.grid.major = element_line(color = "grey99", linewidth = 0.5),
        panel.grid.minor = element_line(color = "white", linewidth = 0.25),
        panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.5))+
  scale_color_manual(values=c("red","black"))



p



final_plot <- plot_grid(
  custom_legend,   # Extracted legend
  p,               # Your main ggplot figure
  ncol = 1,        # Arrange them vertically
  rel_heights = c(0.1,1)  # Give more space to the main plot
)

print(final_plot)


ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Pi_all_autosomes_neoZ_breakpoints.png", 
       plot = final_plot, dpi = 350, width = 14, height = 12, units = "in")


##########################################################################
table(data_filt$in_interval)
data_filt$region <- "Autosomes"
data_filt$region[data_filt$in_interval==TRUE]<-"MORs"
data_filt$region[data_filt$Chr=="Neo_Z"]<-"neo-Z"
table(data_filt$region)

data_filt$region <- factor(data_filt$region,levels=c("neo-Z","Autosomes","MORs"))

# remove the slightly zero observations
#data_filt <- data_[data$V5>0,]

p <- ggplot(data_filt,aes(region,Mean_pi,fill=population))+
  geom_boxplot(outliers=FALSE,notch=TRUE)+
  xlab("")+
#  labs(y = expression(italic(F)[ST]))+
 # stat_compare_means(aes(group=population),comparisons = list(c("Autosomes","MORs"),c("MORs","neo-Z"),c("Autosomes","neo-Z")), 
#                     method = "t.test",
#                     label = "p.signif") +
  theme_minimal()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        legend.position="top")
p



############################################################################################################
# Just the chromosomes with outliers

data_filt_outliers <- data_filt %>%
  filter(Chr %in% c("Chr_4","Chr_6","Chr_8","Chr_9","Chr_10","Chr_13","Chr_14","Chr_15","Chr_17","Chr_18","pri#ptg000015l","pri#ptg000045l"))

y_top <- max(data_filt_outliers$Mean_pi, na.rm = TRUE) * 1.01

p <- ggplot(data_filt_outliers,aes(start/1000000,Mean_pi,color=population))+
  geom_point(alpha=0.35)+
  geom_segment(data=intervals_df,aes(x=start/1e6,xend=end/1e6),y=y_top,yend=y_top,color='blue',linewidth=2,inherit.aes=FALSE)+
  xlab("Chromosome position (Mb)")+
  labs(y = expression(italic(theta[pi])))+
  geom_hline(yintercept = mean_pi, linetype = "dashed", color = "red", size = 1)+
  facet_wrap(.~Chr,scales="free_x")+
  theme_minimal()+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="none",
        legend.title=element_blank(),
        legend.text=element_text(size=16),
          panel.grid.major = element_line(color = "grey99", linewidth = 0.5),
        panel.grid.minor = element_line(color = "white", linewidth = 0.25),
        panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.5))+
  scale_color_manual(values=c("red","black"))


p



final_plot <- plot_grid(
  custom_legend,   # Extracted legend
  p,               # Your main ggplot figure
  ncol = 1,        # Arrange them vertically
  rel_heights = c(0.1,1)  # Give more space to the main plot
)

print(final_plot)


ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Pi_outliers_breakpoints.png", 
       plot = final_plot, dpi = 300, width = 12, height = 8, units = "in")


########################################################################
# box plot of each outlier chromosome

data_filt_outliers$Chrom <- as.character(data_filt_outliers$Chr)
data_filt_outliers$Chrom2 <- as.character(data_filt_outliers$Chr)

data_filt_outliers$Chrom[data_filt_outliers$Chrom2=="Chr_8" & data_filt_outliers$start<10724985] <- "Chr_8.2"
data_filt_outliers$Chrom[data_filt_outliers$Chrom2=="Chr_8" & data_filt_outliers$start>37874992 & data_filt_outliers$start<55974910] <- "Chr_8.1"

table(data_filt_outliers$Chrom)
table(data_filt_outliers$Chr)

data_filt_outliers$Chrom <- factor(data_filt_outliers$Chrom,levels=c("Neo_Z","Chr_1","Chr_3","Chr_4","Chr_6","Chr_7","Chr_8.1","Chr_8.2",
                                     "Chr_9","Chr_10","Chr_11","Chr_12","Chr_13",
                                     "Chr_14","Chr_15","Chr_16","Chr_17","Chr_18",
                                     "Chr_19","Chr_20","Chr_21","Chr_22","Chr_23",
                                     "Chr_24","Chr_25","Chr_26","Chr_27","Chr_28",
                                     "Chr_29","Chr_30","Chr_31","Chr_32","Chr_33",
                                     "Chr_34","Chr_35","Chr_36","Chr_37","Chr_38",
                                     "Chr_39","pri#ptg000015l","pri#ptg000031l",
                                     "pri#ptg000045l","pri#ptg000050l",
                                     "pri#ptg000051l","pri#ptg000085l","pri#ptg000090l"))



data_filt_outliers$popf <- as.factor(data_filt_outliers$population)

p <- ggplot(data_filt_outliers[data_filt_outliers$in_interval==TRUE,],aes(popf,Mean_pi,fill=popf))+
  geom_boxplot(notch=FALSE,outliers=TRUE)+
  facet_wrap(.~Chrom,scales="free")+
  labs(y = expression(italic(theta[pi])))+
  xlab('')+
  stat_compare_means(comparisons=list(c("arid","temperate")),
                     method = "t.test",
                     label = "p.signif",
                     label.y=0.01,
                     tip.length=0,
                     bracket.size=0,
                     color='black',
                     hide.ns=FALSE)+
  theme_minimal()+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="none",
        legend.title=element_blank(),
        legend.text=element_text(size=16))+
  scale_fill_manual(values=c("red","black"))


p

ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Pi_boxplot_outliers.png", 
       plot = p, dpi = 300, width = 12, height = 8, units = "in")



########################################################################
### Look at Tajima's D in outliers

y_top <- max(data_filt_outliers$Tajima, na.rm = TRUE) * 1.01

p <- ggplot(data_filt_outliers,aes(start/1000000,Tajima,color=population))+
  geom_point(alpha=0.35)+
  geom_segment(data=intervals_df,aes(x=start/1e6,xend=end/1e6),y=y_top,yend=y_top,color='blue',linewidth=2,inherit.aes=FALSE)+
  xlab("Chromosome position (Mb)")+
  labs(y = expression(paste("Tajima's ", italic(D))))+
  #geom_hline(yintercept = mean_pi, linetype = "dashed", color = "red", size = 1)+
  facet_wrap(.~Chr,scales="free_x")+
  theme_minimal()+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="none",
        legend.title=element_blank(),
        legend.text=element_text(size=16),
        panel.grid.major = element_line(color = "grey99", linewidth = 0.5),
        panel.grid.minor = element_line(color = "white", linewidth = 0.25),
        panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.5))+
  scale_color_manual(values=c("red","black"))

p




final_plot <- plot_grid(
  custom_legend,   # Extracted legend
  p,               # Your main ggplot figure
  ncol = 1,        # Arrange them vertically
  rel_heights = c(0.1,1)  # Give more space to the main plot
)

print(final_plot)



ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Tajima_outliers_breakpoints.png", 
       plot = final_plot, dpi = 300, width = 12, height = 8, units = "in")


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

data_filt$Interval[data_filt$in_interval==TRUE]<-"MORs"
data_filt$Interval[data_filt$in_interval==FALSE]<-"Autosomes"
data_filt$Interval[data_filt$Chr=="Neo_Z"] <- "neo-Z"

data_filt$Intervalf <- factor(data_filt$Interval,levels=c("neo-Z","Autosomes","MORs"))


data_filt <- data_filt[!is.na(data_filt$population),]
data_filt <- data_filt[!is.na(data_filt$in_interval),]


data_filt$populationf <- droplevels(data_filt$populationf)
data_filt$Intervalf   <- droplevels(data_filt$Intervalf)

p <- ggplot(data_filt,aes(Intervalf,Tajima,fill=population))+
  stat_compare_means(aes(group=population), 
                     method = "t.test",
                     label = "p.signif")+ 
  geom_boxplot(position = position_dodge(width = 0.8),outliers=FALSE)+
  xlab("")+
  labs(y = expression(paste("Tajima's ", italic(D))))+
  theme_minimal()+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="bottom",
        legend.title=element_blank(),
        legend.text=element_text(size=16))+
  scale_fill_manual(values=c("red","black"))



# run the t.test manually
#t.test(V5 ~ region, data = data[data$region %in% c("Autosomes", "neo-Z"), ])
#t.test(V5 ~ region, data = data[data$region %in% c("Autosomes", "MORs"), ])
#t.test(V5 ~ region, data = data[data$region %in% c("neo-Z", "MORs"), ])



p

ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Tajima_autosomes_neoZ_boxplots.png", 
       plot = p, dpi = 300, width = 6, height = 6, units = "in")


##################################################
# Boxplot of Tajima's D in MORs specifically

p <- ggplot(data_filt[data_filt$in_interval==TRUE,],aes(population,Tajima))+
  stat_compare_means(aes(group=population), 
                     method = "t.test",
                     label = "p.signif")+ 
  facet_wrap(.~Chr,scales="free")+
  geom_boxplot(position = position_dodge(width = 0.8),outliers=FALSE)+
  xlab("")+
  labs(y = expression(paste("Tajima's ", italic(D))))+
  theme_minimal()+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="bottom",
        legend.title=element_blank(),
        legend.text=element_text(size=16))+
  scale_fill_manual(values=c("red","black"))


p


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

# Boxplot in temperate pop only 

data_filt$Category=NULL
#data_filt$Category="Autosomes"
data_filt$Chrom <- as.character(data_filt$Chr)
data_filt$Region <- as.character(data_filt$region)
data_filt$Category[data_filt$Region=="MORs" & data_filt$Chrom=="pri#ptg000015l"] = "pri#ptg000015l"
data_filt$Category[data_filt$Region=="MORs" & data_filt$Chrom=="pri#ptg000045l"] = "pri#ptg000045l"
data_filt$Category[data_filt$Region=="MORs" & data_filt$Chrom=="Chr_18"] = "Chr_18"
data_filt$Category[data_filt$Region=="MORs" & data_filt$Chrom=="Chr_17"] = "Chr_17"
data_filt$Category[data_filt$Region=="MORs" & data_filt$Chrom=="Chr_15"] = "Chr_15"
data_filt$Category[data_filt$Region=="MORs" & data_filt$Chrom=="Chr_14"] = "Chr_14"
data_filt$Category[data_filt$Region=="MORs" & data_filt$Chrom=="Chr_13"] = "Chr_13"
data_filt$Category[data_filt$Region=="MORs" & data_filt$Chrom=="Chr_10"] = "Chr_10"
data_filt$Category[data_filt$Region=="MORs" & data_filt$Chrom=="Chr_9"] = "Chr_9"
data_filt$Category[data_filt$Region=="MORs" & data_filt$Chrom=="Chr_8"] = "Chr_8"
data_filt$Category[data_filt$Chrom=="Chr_8" & data_filt$Region=="MORs" & data_filt$start<10724985] <- "Chr_8.2"
data_filt$Category[data_filt$Chrom=="Chr_8" & data_filt$Region=="MORs" & data_filt$start>37874992 & data_filt$start<55974910] <- "Chr_8.1"


data_filt$Category[data_filt$Region=="MORs" & data_filt$Chrom=="Chr_6"] = "Chr_6"
data_filt$Category[data_filt$Region=="MORs" & data_filt$Chrom=="Chr_4"] = "Chr_4"
data_filt$Category[data_filt$Region=="Autosomes"] = "Autosomes"
data_filt$Category[data_filt$Region=="neo-Z"] = "neo-Z"
table(data_filt$Region)
table(data_filt$Category)


str(data_filt)


data_filt_temperate <- data_filt[data_filt$population=="temperate",]
mean(data_filt_temperate$Tajima)
mean(data_filt_temperate$Tajima[data_filt_temperate$Region=="Autosomes"])
median(data_filt_temperate$Tajima[data_filt_temperate$Region=="Autosomes"])
# Create comparison variable
# Define your reference category
my_ref <- "Autosomes"

# Identify all categories except the reference
all_cats <- unique(data_filt_temperate$Category)
compare_cats <- setdiff(all_cats, my_ref)

# For each category in compare_cats, keep only that category and the reference
# and label it as, e.g., "RefCat vs CatB", etc.

library(dplyr)
library(purrr)
library(ggplot2)
df_compare <- map_dfr(compare_cats, function(cat_x) {
  tmp <- data_filt_temperate %>% 
    filter(Category %in% c(my_ref, cat_x)) %>%
    mutate(Comparison = paste(my_ref, "vs", cat_x))
  tmp
})


### Autosome category replicated
p <- ggplot(df_compare,aes(Category,Tajima))+
  stat_compare_means(aes(group=Category),comparisons = list(c("Autosomes","Chr_15"),c("Autosomes","Chr_17"),c("Autosomes","neo-Z"),
                                                            c("Autosomes","Chr_13"),c("Autosomes","Chr_6"),c("Autosomes","Chr_9")), 
                     method = "t.test",
                     label = "p.signif") +
  geom_boxplot(position = position_dodge(width = 0.8),outliers=TRUE,notch=TRUE)+
  geom_hline(yintercept = -1.348771, linetype = "dashed", color = "red", size = 1)+
  facet_wrap(.~Comparison,scales="free")+
  xlab("")+
  labs(y = expression(italic(Tajima)))+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="bottom")
#scale_color_manual(values=c("black","red"))+


p


data_filt_temperate$Category <- factor(data_filt_temperate$Category,levels=c("Autosomes","Chr_4","Chr_6","Chr_8.1","Chr_8.2","Chr_9",
                                                                             "Chr_10","Chr_13","Chr_14","Chr_15","Chr_17",
                                                                             "Chr_18","pri#ptg000015l",
                                                                             "pri#ptg000045l","neo-Z"))

### Single Autosome category
p <- ggplot(data_filt_temperate,aes(Category,Tajima))+
  stat_compare_means(aes(group=Category),comparisons = list(c("Autosomes","Chr_15"),c("Autosomes","Chr_17"),c("Autosomes","neo-Z"),
                                                            c("Autosomes","Chr_13"),c("Autosomes","Chr_6"),c("Autosomes","Chr_9"),
                                                            c("Autosomes","pri#ptg000045l")), 
                       method = "t.test",
                       label = "p.signif") +
  geom_boxplot(position = position_dodge(width = 0.8),outliers=TRUE,notch=TRUE)+
  geom_hline(yintercept = -1.348771, linetype = "dashed", color = "red", size = 1)+
  xlab("")+
  labs(y = expression(italic(Tajima)))+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="bottom")
#scale_color_manual(values=c("black","red"))+


p



### Single Autosome category, significantly lower
p <- ggplot(data_filt_temperate,aes(Category,Tajima,fill=Region))+
  stat_compare_means(aes(group=Category),comparisons = list(c("Autosomes","Chr_15"),
                                                            c("Autosomes","Chr_13"),c("Autosomes","Chr_9"),c("Autosomes","Chr_6")), 
                     method = "t.test",
                     label = "p.signif",
                     label.y= c(-0.5,-0.25,0,0.25,0.5)) +
  geom_boxplot(position = position_dodge(width = 0.8),outliers=TRUE,notch=TRUE)+
  geom_hline(yintercept = -1.325, linetype = "dotted", color = "red", size = 0.5)+
  xlab("")+
  labs(y = expression(paste("Tajima's ", italic(D))))+
  theme_minimal()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=13),
        strip.text = element_text(size=13),
        axis.text.x=element_text(angle=45),
        legend.position="top",
        legend.title=element_blank(),
        legend.text=element_text(size=14))+
scale_fill_manual(values=c("black","blue","white"))


p


ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Tajima_temperate_autos_mors_boxplots.png", 
       plot = p, dpi = 350, width = 12, height = 8, units = "in")





# run the t.test manually
t.test(Tajima ~ Category, data = data_filt_temperate[data_filt_temperate$Category %in% c("Autosomes", "Chr_15"), ])
t.test(Tajima ~ Category, data = data_filt_temperate[data_filt_temperate$Category %in% c("Autosomes", "Chr_6"), ])
t.test(Tajima ~ Category, data = data_filt_temperate[data_filt_temperate$Category %in% c("Autosomes", "Chr_9"), ])











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
                     method = "wilcox.test",
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




