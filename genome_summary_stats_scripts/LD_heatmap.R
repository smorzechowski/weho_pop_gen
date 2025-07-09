# Make LD heatmaps from ngsLD

library(ggplot2)
library(reshape2)
library(stringr)

setwd("~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/ngsLD")
# save as PNGs 500x400


###########################################################################
#data <- read.table('windows_sorted_reformatted_priptg000045l_lea_thin1kb_100kb.ld',header=T,comment.char = '')
#head(data)
#ld_matrix <- acast(data, win1 ~ win2, value.var="mean_r2")
#ld_melt <- melt(ld_matrix, na.rm=TRUE)

#ggplot(ld_melt, aes(x=Var1, y=Var2, fill=value)) +
#  geom_tile() +
#  scale_fill_gradient(low="lightgreen", high="red") +
#  labs(x="Position", y="Position", fill="LD (r²)") +
#  theme_minimal()

###########################################################################
#data <- read.table('windows_sorted_reformatted_priptg000015l_lea_thin1kb_100kb.ld',header=T,comment.char = '')
#ld_matrix <- acast(data, win1 ~ win2, value.var="mean_r2")
#ld_melt <- melt(ld_matrix, na.rm=TRUE)

#ggplot(ld_melt[ld_melt$value<0.2,], aes(x=Var1, y=value)) +
#  geom_point()

###########################################################################
datap15 <- read.table('priptg000015l_lea_thin1kb_final_windows.ld',header=T,comment.char = '')

ld_matrix <- acast(datap15, win1 ~ win2, value.var="mean_r2")
ld_melt <- melt(ld_matrix, na.rm=TRUE)

p <- ggplot(ld_melt, aes(x=Var1/1000000, y=Var2/1000000, fill=value)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="#3F0071",limits=c(0,0.1),oob = scales::squish) +
  labs(x="Chromosome Position (Mb)", y="Chromosome Position (Mb)", fill="LD (r²)") +
  ggtitle("Contig priptg000015l")+
  theme_minimal()+
  theme(axis.title=element_text(size=15),
        axis.text=element_text(size=14),
        legend.text = element_text(size=14),
        legend.title=element_text(size=15))

p 

# Save the plot as a PNG with 300 dpi, for example 10 x 7 inches:
ggsave("Contig_ptg000015l_LD.svg", plot = p, dpi = 300, width = 6, height = 4, units = "in")

###########################################################################
datap45 <- read.table('priptg000045l_lea_thin1kb_final_windows.ld',header=T,comment.char = '')

ld_matrix <- acast(datap45, win1 ~ win2, value.var="mean_r2")
ld_melt <- melt(ld_matrix, na.rm=TRUE)

p <- ggplot(ld_melt, aes(x=Var1/1000000, y=Var2/1000000, fill=value)) +
  geom_tile() +
 # scale_fill_gradient(low="white", high="#3F0071") +
 scale_fill_gradient(low="white", high="#3F0071",limits=c(0,0.15),oob = scales::squish) +
  labs(x="Chromosome Position (Mb)", y="Chromosome Position (Mb)", fill="LD (r²)") +
  ggtitle("Contig priptg000045l")+
  theme_minimal()+
  theme(axis.title=element_text(size=15),
        axis.text=element_text(size=14),
        legend.text = element_text(size=14),
        legend.title=element_text(size=15))

p 

# Save the plot as a PNG with 300 dpi, for example 10 x 7 inches:
ggsave("Contig_ptg000045l_LD.svg", plot = p, dpi = 300, width = 6, height = 4, units = "in")

###########################################################################
data4 <- read.table('scaffold_4_RagTag_lea_thin1kb_final_windows.ld',header=T,comment.char = '')

ld_matrix <- acast(data4, win1 ~ win2, value.var="mean_r2")
ld_melt <- melt(ld_matrix, na.rm=TRUE)

p <- ggplot(ld_melt, aes(x=Var1/1000000, y=Var2/1000000, fill=value)) +
  geom_tile() +
  #scale_fill_gradient(low="white", high="#3F0071") +
  scale_fill_gradient(low="white", high="#3F0071",limits=c(0,0.15),oob = scales::squish) +
  labs(x="Chromosome Position (Mb)", y="Chromosome Position (Mb)", fill="LD (r²)") +
  ggtitle("Chromosome 4")+
  theme_minimal()+
  theme(axis.title=element_text(size=15),
        axis.text=element_text(size=14),
        legend.text = element_text(size=14),
        legend.title=element_text(size=15))

# Save the plot as a PNG with 300 dpi, for example 10 x 7 inches:
ggsave("Chr_4_LD.svg", plot = p, dpi = 300, width = 6, height = 4, units = "in")




###########################################################################
data6 <- read.table('scaffold_6_RagTag_lea_thin1kb_final_windows.ld',header=T,comment.char = '')

ld_matrix <- acast(data6, win1 ~ win2, value.var="mean_r2")
ld_melt <- melt(ld_matrix, na.rm=TRUE)

p <- ggplot(ld_melt, aes(x=Var1/1000000, y=Var2/1000000, fill=value)) +
  geom_tile() +
  #scale_fill_gradient(low="white", high="#3F0071") +
  scale_fill_gradient(low="white", high="#3F0071",limits=c(0,0.1),oob = scales::squish) +
  labs(x="Chromosome Position (Mb)", y="Chromosome Position (Mb)", fill="LD (r²)") +
  ggtitle("Chromosome 6")+
  theme_minimal()+
  theme(axis.title=element_text(size=15),
        axis.text=element_text(size=14),
        legend.text = element_text(size=14),
        legend.title=element_text(size=15))

p
# Save the plot as a PNG with 300 dpi, for example 10 x 7 inches:
ggsave("Chr_6_LD.svg", plot = p, dpi = 300, width = 6, height = 4, units = "in")

###########################################################################
data8 <- read.table('scaffold_8_RagTag_lea_thin1kb_final_windows.ld',header=T,comment.char = '')

ld_matrix <- acast(data8, win1 ~ win2, value.var="mean_r2")
ld_melt <- melt(ld_matrix, na.rm=TRUE)

p <- ggplot(ld_melt, aes(x=Var1/1000000, y=Var2/1000000, fill=value)) +
  geom_tile() +
  #scale_fill_gradient(low="white", high="#3F0071") +
  scale_fill_gradient(low="white", high="#3F0071",limits=c(0,0.1),oob = scales::squish) +
  labs(x="Chromosome Position (Mb)", y="Chromosome Position (Mb)", fill="LD (r²)") +
  ggtitle("Chromosome 8")+
  theme_minimal()+
  theme(axis.title=element_text(size=15),
        axis.text=element_text(size=14),
        legend.text = element_text(size=14),
        legend.title=element_text(size=15))
p

# Save the plot as a PNG with 300 dpi, for example 10 x 7 inches:
ggsave("Chr_8_LD.svg", plot = p, dpi = 300, width = 6, height = 4, units = "in")

###########################################################################
data9 <- read.table('scaffold_9_RagTag_lea_thin1kb_final_windows.ld',header=T,comment.char = '')

ld_matrix <- acast(data9, win1 ~ win2, value.var="mean_r2")
ld_melt <- melt(ld_matrix, na.rm=TRUE)

p <- ggplot(ld_melt, aes(x=Var1/1000000, y=Var2/1000000, fill=value)) +
  geom_tile() +
  #scale_fill_gradient(low="white", high="#3F0071") +
  scale_fill_gradient(low="white", high="#3F0071",limits=c(0,0.1),oob = scales::squish) +
  labs(x="Chromosome Position (Mb)", y="Chromosome Position (Mb)", fill="LD (r²)") +
  ggtitle("Chromosome 9")+
  theme_minimal()+
  theme(axis.title=element_text(size=15),
        axis.text=element_text(size=14),
        legend.text = element_text(size=14),
        legend.title=element_text(size=15))

# Save the plot as a PNG with 300 dpi, for example 10 x 7 inches:
ggsave("Chr_9_LD.svg", plot = p, dpi = 300, width = 6, height = 4, units = "in")

###########################################################################
data10 <- read.table('scaffold_10_RagTag_lea_thin1kb_final_windows.ld',header=T,comment.char = '')

ld_matrix <- acast(data10, win1 ~ win2, value.var="mean_r2")
ld_melt <- melt(ld_matrix, na.rm=TRUE)

p <- ggplot(ld_melt, aes(x=Var1/1000000, y=Var2/1000000, fill=value)) +
  geom_tile() +
  #scale_fill_gradient(low="white", high="#3F0071") +
  scale_fill_gradient(low="white", high="#3F0071",limits=c(0,0.3),oob = scales::squish) +
  labs(x="Chromosome Position (Mb)", y="Chromosome Position (Mb)", fill="LD (r²)") +
  ggtitle("Chromosome 10")+
  theme_minimal()+
  theme(axis.title=element_text(size=15),
        axis.text=element_text(size=14),
        legend.text = element_text(size=14),
        legend.title=element_text(size=15))

p

# Save the plot as a PNG with 300 dpi, for example 10 x 7 inches:
ggsave("Chr_10_LD.svg", plot = p, dpi = 300, width = 6, height = 4, units = "in")

###########################################################################
data13 <- read.table('scaffold_13_RagTag_lea_thin1kb_final_windows.ld',header=T,comment.char = '')

ld_matrix <- acast(data13, win1 ~ win2, value.var="mean_r2")
ld_melt <- melt(ld_matrix, na.rm=TRUE)

p <- ggplot(ld_melt, aes(x=Var1/1000000, y=Var2/1000000, fill=value)) +
  geom_tile() +
  #scale_fill_gradient(low="white", high="#3F0071") +
  scale_fill_gradient(low="white", high="#3F0071",limits=c(0,0.1),oob = scales::squish) +
  labs(x="Chromosome Position (Mb)", y="Chromosome Position (Mb)", fill="LD (r²)") +
  ggtitle("Chromosome 13")+
  theme_minimal()+
  theme(axis.title=element_text(size=15),
        axis.text=element_text(size=14),
        legend.text = element_text(size=14),
        legend.title=element_text(size=15))

p

# Save the plot as a PNG with 300 dpi, for example 10 x 7 inches:
ggsave("Chr_13_LD.svg", plot = p, dpi = 300, width = 6, height = 4, units = "in")


###########################################################################
data14 <- read.table('scaffold_14_RagTag_lea_thin1kb_final_windows.ld',header=T,comment.char = '')

ld_matrix <- acast(data14, win1 ~ win2, value.var="mean_r2")
ld_melt <- melt(ld_matrix, na.rm=TRUE)

p <- ggplot(ld_melt, aes(x=Var1/1000000, y=Var2/1000000, fill=value)) +
  geom_tile() +
  #scale_fill_gradient(low="white", high="#3F0071") +
  scale_fill_gradient(low="white", high="#3F0071",limits=c(0,0.15),oob = scales::squish) +
  labs(x="Chromosome Position (Mb)", y="Chromosome Position (Mb)", fill="LD (r²)") +
  ggtitle("Chromosome 14")+
  theme_minimal()+
  theme(axis.title=element_text(size=15),
        axis.text=element_text(size=14),
        legend.text = element_text(size=14),
        legend.title=element_text(size=15))


p

# Save the plot as a PNG with 300 dpi, for example 10 x 7 inches:
ggsave("Chr_14_LD.svg", plot = p, dpi = 300, width = 6, height = 4, units = "in")

###########################################################################
data15 <- read.table('scaffold_15_RagTag_lea_thin1kb_final_windows.ld',header=T,comment.char = '')

ld_matrix <- acast(data15, win1 ~ win2, value.var="mean_r2")
ld_melt <- melt(ld_matrix, na.rm=TRUE)

p <- ggplot(ld_melt, aes(x=Var1/1000000, y=Var2/1000000, fill=value)) +
  geom_tile() +
  #scale_fill_gradient(low="white", high="#3F0071") +
  scale_fill_gradient(low="white", high="#3F0071",limits=c(0,0.2),oob = scales::squish) +
  labs(x="Chromosome Position (Mb)", y="Chromosome Position (Mb)", fill="LD (r²)") +
  ggtitle("Chromosome 15")+
  theme_minimal()+
  theme(axis.title=element_text(size=15),
        axis.text=element_text(size=14),
        legend.text = element_text(size=14),
        legend.title=element_text(size=15))

p

# Save the plot as a PNG with 300 dpi, for example 10 x 7 inches:
ggsave("Chr_15_LD.svg", plot = p, dpi = 300, width = 6, height = 4, units = "in")

###########################################################################
data17 <- read.table('scaffold_17_RagTag_lea_thin1kb_final_windows.ld',header=T,comment.char = '')

ld_matrix <- acast(data17, win1 ~ win2, value.var="mean_r2")
ld_melt <- melt(ld_matrix, na.rm=TRUE)

p <- ggplot(ld_melt, aes(x=Var1/1000000, y=Var2/1000000, fill=value)) +
  geom_tile() +
  #scale_fill_gradient(low="lightyellow", high="#3F0071") +
  scale_fill_gradient(low="white", high="#3F0071",limits=c(0,0.1),oob = scales::squish) +
  labs(x="Chromosome Position (Mb)", y="Chromosome Position (Mb)", fill="LD (r²)") +
  ggtitle("Chromosome 17")+
  theme_minimal()+
  theme(axis.title=element_text(size=15),
        axis.text=element_text(size=14),
        legend.text = element_text(size=14),
        legend.title=element_text(size=15))


p

# Save the plot as a PNG with 300 dpi, for example 10 x 7 inches:
ggsave("Chr_17_LD.svg", plot = p, dpi = 300, width = 6, height = 4, units = "in")

###########################################################################
data18 <- read.table('scaffold_18_RagTag_lea_thin1kb_final_windows.ld',header=T,comment.char = '')

ld_matrix <- acast(data18, win1 ~ win2, value.var="mean_r2")
ld_melt <- melt(ld_matrix, na.rm=TRUE)

p <- ggplot(ld_melt, aes(x=Var1/1000000, y=Var2/1000000, fill=value)) +
  geom_tile() +
   #scale_fill_gradient(low="white", high="#3F0071") +
   scale_fill_gradient(low="white", high="#3F0071",limits=c(0,0.2),oob = scales::squish) +
  labs(x="Chromosome Position (Mb)", y="Chromosome Position (Mb)", fill="LD (r²)") +
  theme_minimal()+
  ggtitle("Chromosome 18")+
  theme_minimal()+
  theme(axis.title=element_text(size=15),
        axis.text=element_text(size=14),
        legend.text = element_text(size=14),
        legend.title=element_text(size=15))

p

# Save the plot as a PNG with 300 dpi, for example 10 x 7 inches:
ggsave("Chr_18_LD.svg", plot = p, dpi = 300, width = 6, height = 4, units = "in")


########################################################################################
# Combine all data and make a facet grid plot

data_all <- rbind(datap45,datap15,data4,data6,data8,data9,data10,data13,data14,data15,data17,data18)


data_all$chr <- str_replace(data_all$chr,'scaffold','Chr')
data_all$chr <- str_remove(data_all$chr,'_RagTag')

data_all$chr <- factor(data_all$chr,levels=c("Chr_4","Chr_6","Chr_8","Chr_9","Chr_10","Chr_13","Chr_14","Chr_15","Chr_17","Chr_18","pri#ptg000015l","pri#ptg000045l"))




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


intervals_df$Chr <- factor(intervals_df$Chr,levels=c("Chr_4","Chr_6",
                                                     "Chr_8","Chr_9","Chr_10",
                                                     "Chr_13","Chr_14",
                                                     "Chr_15","Chr_17","Chr_18","pri#ptg000015l",
                                                     "pri#ptg000045l"))

# Compute max for each chromosome
max_vals <- tapply(data_all$win2, data_all$chr, max)

# Convert to a data frame with two columns: 'chr' and 'max'
y_top <- data.frame(Chr = names(max_vals), y_top = as.vector(max_vals)*1.1)

intervals_df <- data.frame(left_join(intervals_df,y_top))
colnames(intervals_df)[1] <- "chr"
head(intervals_df)

intervals_df$chr <- factor(intervals_df$chr,levels=c("Chr_4","Chr_6",
                                                     "Chr_8","Chr_9","Chr_10",
                                                     "Chr_13","Chr_14",
                                                     "Chr_15","Chr_17","Chr_18","pri#ptg000015l",
                                                     "pri#ptg000045l"))



p <- ggplot(data_all, aes(x=win1/1000000, y=win2/1000000, fill=mean_r2)) +
  geom_tile() +
  #scale_fill_gradient(low="white", high="#3F0071") +
  scale_fill_gradient(low="white", high="#3F0071",limits=c(0,0.09),oob = scales::squish) +
  geom_segment(data=intervals_df,aes(x=start/1e6,xend=end/1e6,y=y_top/1e6,yend=y_top/1e6,group=chr),color='#3F0071',linewidth=2,inherit.aes = FALSE)+
  facet_wrap(.~chr,scales="free")+
  labs(x="Chromosome Position (Mb)", y="Chromosome Position (Mb)", fill="LD (r²)") +
  theme_minimal()+
 # ggtitle("Chromosome 18")+
  theme_minimal()+
  theme(axis.title=element_text(size=15),
        axis.text=element_text(size=14),
        legend.text = element_text(size=14),
        legend.title=element_text(size=15),
        strip.text = element_text(size=14))

p

# Save the plot as a PNG with 300 dpi, for example 10 x 7 inches:
ggsave("all_chromosomes_LD.png", plot = p, dpi = 500, width = 10, height = 7, units = "in")

