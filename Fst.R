# Plot Fst across the chromosomes
# Sophie C. M. Orzechowski
# February 2025


######################################################################################
library(stringr)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(gghalves)
library(cowplot)
library(gridExtra)
library(ggrepel)
setwd("~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/fst")
data <- read.table("temperate_arid.pbs.fst.txt",header=F,comment.char = "",skip=1)
data <- read.table("temperate_arid.pbs.fst_50kb.txt",header=F,comment.char = "",skip=1)

Zdata <- read.table("temperate_arid_neoZ_males.pbs.fst_50kb.txt",header=F,comment.char = "",skip=1)


######################################################################################

data$Chr <- str_replace(data$V2,'scaffold','Chr')
data$Chr <- str_remove(data$Chr,'_RagTag')


Zdata$Chr <- "Neo_Z"

#Zdata$Chr <- str_replace(Zdata$V2,'scaffold','Chr')
#Zdata$Chr <- str_remove(Zdata$Chr,'_RagTag')

data <- rbind(data,Zdata)

data$start <- data$V3 - 25000

# for combine_summary_metrics.R
#colnames(data)[5] <- 'Fst'
#Fstdata <- data





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



y_top <- max(data$V5, na.rm = TRUE) * 1.01


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


######################################################################################
# All chromosomes


# Plot vlines for neo-Z
newPAR_vline <- data.frame(Chr = "Neo_Z",xcoord=111.6)
newPAR_vline$Chr <- as.factor(newPAR_vline$Chr)

# Plot vlines for neo-Z
addedZ_vline <- data.frame(Chr = "Neo_Z",xcoord=74.4)
addedZ_vline$Chr <- as.factor(addedZ_vline$Chr)

# Plot intervals for neo-Z
intervals_neoZ <- data.frame(Chr=c("Neo_Z","Neo_Z","Neo_Z"),start=c(0,74.4,111.6),end=c(74.4,111.6,128.9))
intervals_neoZ$Chr <- as.factor(intervals_neoZ$Chr)

p <- ggplot(data,aes(V3/1000000,V5,color=in_interval))+
  geom_point()+
  geom_segment(data=intervals_df,aes(x=start/1e6,xend=end/1e6),y=y_top,yend=y_top,color='blue',linewidth=2,inherit.aes=FALSE)+
  # Black outline segments (slightly thicker)
  geom_segment(data = intervals_neoZ,aes(x=start,xend=end),y = y_top, yend = y_top,color = 'black', linewidth = 2.8,inherit.aes=FALSE) +
  geom_segment(data=intervals_neoZ,aes(x=start,xend=end),y=y_top,yend=y_top,color=c('black','white','#857E13'),linewidth=2,inherit.aes=FALSE)+
 ylim(0,0.8)+
  xlab("Chromosome position (Mb)")+
  labs(y = expression(italic(F)[ST]))+
  facet_wrap(.~Chr,scales="free_x")+
#  geom_vline(data = newPAR_vline,aes(xintercept = xcoord),color = "red",linewidth = 1,linetype='dotted')+
#  geom_vline(data = addedZ_vline,aes(xintercept = xcoord),color = "black",linewidth = 1,linetype='dotted')+
  theme_minimal()+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="none",
          panel.grid.major = element_line(color = "grey99", linewidth = 0.5),
          panel.grid.minor = element_line(color = "white", linewidth = 0.25),
        panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.5)
        )+
  #scale_color_manual(values=c("darkgrey","darkblue"))+
  scale_color_manual(values=c("#696969","#696969"))


p

final_plot <- plot_grid(
  custom_legend,   # Extracted legend
  p,               # Your main ggplot figure
  ncol = 1,        # Arrange them vertically
 rel_heights = c(0.1,1)  # Give more space to the main plot
)

print(final_plot)


#ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Fst_all_chromosomes_neoZ_breakpoints.png", 
#       plot = final_plot, dpi = 350, width = 14, height = 12, units = "in")

ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Fst_all_chromosomes_neoZ_breakpoints_spelled_out.png", 
       plot = final_plot, dpi = 350, width = 14, height = 12, units = "in")



##################################################################################
# Focus on just the neo-Z and add gene annotations of significant GEA hits

GEA <- read.table("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/GEA_hits_neoZ.txt",sep=" ",head=F)
GEA_order <- GEA[order(GEA$V3),]
GEA_order_nondup <- GEA_order[!duplicated(GEA_order$V3),]

y_top <- max(neoZ$V5, na.rm = TRUE) * 1.1
# Plot intervals for neo-Z
intervals_neoZ <- data.frame(Chr=c("Neo_Z","Neo_Z","Neo_Z"),start=c(0,74.4,111.6),end=c(74.4,111.6,128.9))
intervals_neoZ$Chr <- as.factor(intervals_neoZ$Chr)

neoZ <- data[data$V2=="scaffold_2_RagTag",]

p <- ggplot(neoZ,aes(start/1000000,V5))+
  geom_point()+
  # Black outline segments (slightly thicker)
  geom_segment(data = intervals_neoZ,aes(x=start,xend=end),y = y_top, yend = y_top,color = 'black', linewidth = 5,inherit.aes=FALSE) +
  geom_segment(data=intervals_neoZ,aes(x=start,xend=end),y=y_top,yend=y_top,color=c('black','white','#857E13'),linewidth=4.2,inherit.aes=FALSE)+
  ylim(-0.01,0.25)+
  xlab("Chromosome position (Mb)")+
  labs(y = expression(italic(F)[ST]))+
  theme_minimal()+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="none",
        panel.grid.major = element_line(color = "grey99", linewidth = 0.5),
        panel.grid.minor = element_line(color = "white", linewidth = 0.25),
        panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.5)
  )+
  #scale_color_manual(values=c("darkgrey","darkblue"))+
  scale_color_manual(values=c("#696969","#696969"))+
  geom_text_repel(
    data = GEA,       # Separate annotation data frame
    aes(x = V2/1000000, y = 0.1, label = V3),
    color = "blue",        # Customize text color
    vjust = -1,            # Adjust vertical position
    size = 5,
    max.overlaps = 20
  )


p











######################################################################################
# box plot of Fst inside and outside inversions and on the neo-Z

table(data$in_interval)
data$region <- "Autosomes"
data$region[data$in_interval==TRUE]<-"MORs"
data$region[data$V2=="scaffold_2_RagTag"]<-"neo-Z"
table(data$region)

data$region <- factor(data$region,levels=c("neo-Z","Autosomes","MORs"))

# remove the slightly zero observations
data_filt <- data[data$V5>0,]

p <- ggplot(data_filt,aes(region,V5,group=region))+
  geom_boxplot(outliers=TRUE)+
  #ylim(0,1)+
  xlab("")+
  labs(y = expression(italic(F)[ST]))+
  stat_compare_means(aes(group=region),comparisons = list(c("Autosomes","MORs"),c("MORs","neo-Z"),c("Autosomes","neo-Z")), 
                     method = "t.test",
                     label = "p.signif") +
  theme_minimal()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        legend.position="none")
p


ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Fst_autosomes_neoZ_boxplots.png", 
       plot = p, dpi = 300, width = 6, height = 4, units = "in")



# run the t.test manually
t.test(V5 ~ region, data = data_filt[data_filt$region %in% c("Autosome-wide", "neo-Z"), ])
t.test(V5 ~ region, data = data_filt[data_filt$region %in% c("Autosome-wide", "PCRs"), ])
t.test(V5 ~ region, data = data_filt[data_filt$region %in% c("neo-Z", "PCRs"), ])

# look at ggplot data
ggplot_build(p)$data


data_filt <- data[data$region!="PCRs",]
data_filt$region <- droplevels(data_filt$region)
data_filt <- data_filt[data_filt$V5>0,]

######################################################################################
# remove the PCRs to look at the difference between Autosome-wide and neo-Z
p <- ggplot(data_filt,aes(region,V5))+
  geom_boxplot(outliers=FALSE,notch=TRUE)+
 # ylim(0,1)+
  xlab("")+
  labs(y = expression(italic(F)[ST]))+
 # stat_compare_means(comparisons = list(c("Autosome-wide", "neo-Z")), 
#                     method = "t.test",
#                     label = "p.signif") +
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        legend.position="none")
#  scale_x_discrete(labels = c("TRUE" = "PCRs", "FALSE" = "Autosome-wide"))

p
t.test(V5 ~ region, data = data_filt[data_filt$region %in% c("Autosome-wide", "neo-Z"), ])
wilcox.test(V5 ~ region, data = data_filt[data_filt$region %in% c("Autosome-wide", "neo-Z"), ],exact=FALSE)

median(data$V5[data$in_interval==TRUE])
median(data$V5[data$in_interval==FALSE])

median(data_filt$V5[data_filt$in_interval==TRUE])
median(data_filt$V5[data_filt$in_interval==FALSE])
median(data_filt$V5[data_filt$V2=="scaffold_2_RagTag"])

######################################################################################
# Histogram of Fst inside and outside inversions
ggplot(data,aes(V5,fill=in_interval,color=in_interval))+
  geom_histogram(bins=100,alpha=.3,position='identity')+
  ylab("Count (50 kb windows)")+
  labs(x = expression(italic(F)[ST]))+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        legend.position="none")



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

######################################################################################
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

p
# Save the plot as a PNG with 300 dpi, for example 10 x 7 inches:
ggsave("Fst_contig_00047l_breakpoints.svg", plot = p, dpi = 300, width = 10, height = 7, units = "in")



######################################################################################
data_subset <- data[data$V2=="pri#ptg000045l",]
intervals_df_subset <- intervals_df[intervals_df$Chr=="pri#ptg000045l",]

# Individual plots
p <- ggplot(data_subset,aes(V3/1000000,V5))+
  geom_point()+
  geom_line(linewidth=1)+
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
        legend.position="none")+
  scale_x_continuous(breaks = seq(0, 18, by = 2))
p

ggsave("Fst_priptg000045l_breakpoints_updated_inset.png", plot = p, dpi = 400, width = 7, height = 5, units = "in")

# Save the plot as a PNG with 300 dpi, for example 10 x 7 inches:
#ggsave("Fst_priptg000015l_breakpoints.svg", plot = p, dpi = 300, width = 10, height = 7, units = "in")


######################################################################################
data_subset <- data[data$V2=="pri#ptg000015l",]
intervals_df_subset <- intervals_df[intervals_df$Chr=="pri#ptg000015l",]

# Individual plots
p <- ggplot(data_subset,aes(V3/1000000,V5))+
  geom_point()+
  geom_line(linewidth=1)+
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
        legend.position="none")+
  scale_x_continuous(breaks = seq(0, 18, by = 2))
p

ggsave("Fst_priptg000015l_breakpoints_updated_inset.png", plot = p, dpi = 400, width = 7, height = 5, units = "in")

# Save the plot as a PNG with 300 dpi, for example 10 x 7 inches:
#ggsave("Fst_priptg000015l_breakpoints.svg", plot = p, dpi = 300, width = 10, height = 7, units = "in")



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
  geom_vline(xintercept = 2.8, color = "black", linetype = "dashed") +
  geom_vline(xintercept = 22.9, color = "black", linetype = "dashed") +
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
ggsave("Fst_scaffold_13_breakpoints_colors_updated_inset.svg", plot = p, dpi = 350, width = 7, height = 5, units = "in")
#ggsave("Fst_scaffold_13_breakpoints.svg", plot = p, dpi = 350, width = 10, height = 7, units = "in")


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




##########################################################################################################
# Fst on the neo-Z

data <- read.table("temperate_arid_neoZ_males.pbs.fst_50kb.txt",header=F,comment.char = "",skip=1)

y_top <- max(data$V5, na.rm = TRUE) * 1.1

data <- data[data$V5>0,]

p <- ggplot(data,aes(V3/1000000,V5))+
  geom_point()+
#  ylim(0,0.8)+
  xlab("Chromosome position (Mb)")+
  labs(y = expression(italic(F)[ST]))+
  geom_vline(xintercept = 74.4,
    color = "darkgrey",
    linewidth = 1,
    linetype ='dashed')+
  labs(y = expression(italic(F)[ST]))+
  geom_vline(xintercept = 111.6,
             color = "red",
             linewidth = 1,
             linetype ='dashed')+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="none")


p

# Fst estimates

# ancestral Z
mean(data$V5[data$V3<74400000])
median(data$V5[data$V3<74400000])

# added Z
mean(data$V5[data$V3>74400000 & data$V3 < 111600000])
median(data$V5[data$V3>74400000 & data$V3 < 111600000])

# new PAR
mean(data$V5[data$V3 > 111600000])
median(data$V5[data$V3 > 111600000])

ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Fst_neoZ_regions_lines.png", 
       plot = p, dpi = 300, width = 6, height = 6, units = "in")

