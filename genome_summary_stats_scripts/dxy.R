# calculate Dxy
# Sophie M Orzechowski
# March 2025

#######################################################################################

library(ggplot2)
library(stringr)
library(dplyr)

setwd("~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/dxy")


#data <- read.table('dxy_50kb_alldata.txt',comment.char = "",quote="",header=T)
data <- read.table('dxy_50kb_alldata_fixed.txt',comment.char = "",quote="",header=T)
#Zdata <- read.table('dxy_50kb_scaffold_2_neoZ.txt',header=T,quote="",comment.char = "")

#data <- rbind(data,Zdata)
mean(data$dxy)



Dxydata <- data



data$Chr <- str_replace(data$scaffold,'scaffold','Chr')
data$Chr <- str_remove(data$Chr,'_RagTag')

data$Chr[data$scaffold=="scaffold_2_RagTag"]<-"Neo_Z"
table(data$Chr)

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


###########################################################################

data_filt <- data[data$dxy<0.08921566,]
data_filt <- data[data$dxy<0.06,]

y_top <- max(data_filt$dxy, na.rm = TRUE) * 1.01



# Plot intervals for neo-Z
intervals_neoZ <- data.frame(Chr=c("Neo_Z","Neo_Z","Neo_Z"),start=c(0,74.4,111.6),end=c(74.4,111.6,128.9))
intervals_neoZ$Chr <- as.factor(intervals_neoZ$Chr)


mean_dxy <- mean(data_filt$dxy)

p <- ggplot(data_filt,aes(start/1000000,dxy,color=in_interval))+
  geom_point()+
  geom_segment(data=intervals_df,aes(x=start/1e6,xend=end/1e6),y=y_top,yend=y_top,color='blue',linewidth=2,inherit.aes=FALSE)+
  # Black outline segments (slightly thicker)
  geom_segment(data = intervals_neoZ,aes(x=start,xend=end),y = y_top, yend = y_top,color = 'black', linewidth = 2.8,inherit.aes=FALSE) +
  geom_segment(data=intervals_neoZ,aes(x=start,xend=end),y=y_top,yend=y_top,color=c('black','white','#857E13'),linewidth=2,inherit.aes=FALSE)+
  xlab("Chromosome position (Mb)")+
  labs(y = expression(italic(D)[xy]))+
  geom_hline(yintercept = mean_dxy, linetype = "dashed", color = "red", size = 1)+
  facet_wrap(.~Chr,scales="free_x")+
  theme_minimal()+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=12),
        strip.text = element_text(size=13),
        legend.position="none",
        panel.grid.major = element_line(color = "grey99", linewidth = 0.5),
        panel.grid.minor = element_line(color = "white", linewidth = 0.25),
        panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.5))+
 # scale_color_manual(values=c("grey","grey"))+
  scale_color_manual(values=c("#696969","#696969"))



p


final_plot <- plot_grid(
  custom_legend,   # Extracted legend
  p,               # Your main ggplot figure
  ncol = 1,        # Arrange them vertically
  rel_heights = c(0.1,1)  # Give more space to the main plot
)

print(final_plot)


ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Dxy_all_chromosomes_neoZ_breakpoints_fixed.png", 
       plot = final_plot, dpi = 350, width = 14, height = 12, units = "in")





# Individual plots
#####################################################################################

#data <- read.table('dxy_50kb_scaffold_6.txt',header=T)
data <- read.table('dxy_50kb_scaffold_6_v3.txt',header=T)

#####################################################################################
xstart <- 8125001
xend <- 32874990
y_top <- 0.003

ggplot(data,aes(start/1000000,dxy))+
  geom_point()+
  ylim(0,.004)+
  stat_smooth(method='loess')+
  geom_segment(x=xstart/1e6,xend=xend/1e6,y=y_top,yend=y_top,color='blue',linewidth=2)



#####################################################################################
data <- read.table('dxy_50kb_scaffold_9_RagTag.txt',header=T)

xstart<-1631342
xend <-18339159
y_top<- 0.004


data_filt <- data[data$dxy<0.003,]

ggplot(data_filt,aes(start/1000000,dxy))+
  geom_point()+
  ylim(0,0.005)+
  stat_smooth(method='loess')+
  geom_segment(x=xstart/1e6,xend=xend/1e6,y=y_top,yend=y_top,color='blue',linewidth=2)



#####################################################################################
data <- read.table('dxy_50kb_scaffold_10.txt',header=T)

xstart<- 9675047
xend <- 22345924
y_top<- 0.004


ggplot(data,aes(start/1000000,dxy))+
  geom_point()+
  ylim(0,0.005)+
  stat_smooth(method='loess')+
  geom_segment(x=xstart/1e6,xend=xend/1e6,y=y_top,yend=y_top,color='blue',linewidth=2)


#####################################################################################
data <- read.table('dxy_50kb_scaffold_11.txt',header=T)

xstart<- NA
xend <- NA
y_top<- NA


ggplot(data,aes(start/1000000,dxy))+
  geom_point()+
 # ylim(0,0.005)+
  stat_smooth(method='loess')
  #geom_segment(x=xstart/1e6,xend=xend/1e6,y=y_top,yend=y_top,color='blue',linewidth=2)



#####################################################################################
data <- read.table('dxy_50kb_scaffold_12.txt',header=T)

xstart<- NA
xend <- NA
y_top<- NA


ggplot(data,aes(start/1000000,dxy))+
  geom_point()+
  ylim(0,0.005)+
  stat_smooth(method='loess')


#####################################################################################
data <- read.table('dxy_50kb_scaffold_13.txt',header=T)

xstart<-15000000
xend <- 23500000
y_top<-0.003

data_filt <- data[data$dxy<0.003,]

ggplot(data_filt,aes(start/1000000,dxy))+
  geom_point()+
  ylim(0,0.004)+
  stat_smooth(method='loess')+
  geom_segment(x=xstart/1e6,xend=xend/1e6,y=y_top,yend=y_top,color='blue',linewidth=2)





#####################################################################################
data <- read.table('dxy_50kb_scaffold_14.txt',header=T)

xstart<- 1
xend <- 10174990
y_top<-0.003


ggplot(data,aes(start/1000000,dxy))+
  geom_point()+
  ylim(0,0.004)+
  stat_smooth(method='loess')+
  geom_segment(x=xstart/1e6,xend=xend/1e6,y=y_top,yend=y_top,color='blue',linewidth=2)




#####################################################################################
data <- read.table('dxy_50kb_scaffold_15.txt',header=T)

xstart<-  1
xend <- 8774970
y_top<-0.003

data_filt <- data[data$dxy<0.003,]

ggplot(data_filt,aes(start/1000000,dxy))+
  geom_point()+
  ylim(0,0.004)+
  stat_smooth(method='loess')+
  geom_segment(x=xstart/1e6,xend=xend/1e6,y=y_top,yend=y_top,color='blue',linewidth=2)


#####################################################################################
data <- read.table('dxy_50kb_scaffold_16.txt',header=T)

xstart<- NA
xend <- NA
y_top<-NA


ggplot(data[data$dxy<0.003,],aes(start/1000000,dxy))+
  geom_point()+
  ylim(0,0.004)+
  stat_smooth(method='loess')
# geom_segment(x=xstart/1e6,xend=xend/1e6,y=y_top,yend=y_top,color='blue',linewidth=2)



#####################################################################################
data <- read.table('dxy_50kb_scaffold_17.txt',header=T)

xstart<-12708285
xend <-19617658
y_top<-0.003
data_filt <- data[data$dxy<0.003,]

ggplot(data_filt,aes(start/1000000,dxy))+
  geom_point()+
  ylim(0,0.004)+
  stat_smooth(method='loess')+
  geom_segment(x=xstart/1e6,xend=xend/1e6,y=y_top,yend=y_top,color='blue',linewidth=2)


#####################################################################################
data <- read.table('dxy_50kb_scaffold_18.txt',header=T)

xstart<- 2865069
xend <- 10487981
y_top<-0.003
data_filt <- data[data$dxy<0.003,]

ggplot(data_filt,aes(start/1000000,dxy))+
  geom_point()+
  ylim(0,0.004)+
  stat_smooth(method='loess')+
  geom_segment(x=xstart/1e6,xend=xend/1e6,y=y_top,yend=y_top,color='blue',linewidth=2)



#####################################################################################
data <- read.table('dxy_50kb_pri#ptg000015l.txt',header=T,quote="",comment.char = "")

xstart<- 6864979
xend <- 17909547
y_top<-0.003
data_filt <- data[data$dxy<0.003,]

ggplot(data_filt,aes(start/1000000,dxy))+
  geom_point()+
  ylim(0,0.004)+
  stat_smooth(method='loess')+
  geom_segment(x=xstart/1e6,xend=xend/1e6,y=y_top,yend=y_top,color='blue',linewidth=2)





#####################################################################################
data <- read.table('dxy_50kb_scaffold_2_neoZ.txt',header=T,quote="",comment.char = "")

xstart<-1
xend <-13100247
y_top<-0.003
data_filt <- data[data$dxy<0.003,]

ggplot(data_filt,aes(start/1000000,dxy))+
  geom_point()+
  #ylim(0,0.004)+
  stat_smooth(method='loess')





