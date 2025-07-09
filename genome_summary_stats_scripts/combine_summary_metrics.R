# Combine summary statistics for Fst, Dxy and theta



library(dplyr)
library(patchwork)
library(cowplot)



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






##################################################################
Chr4Fst <- Fstdata[Fstdata$V2=="scaffold_4_RagTag"&Fstdata$V4>16000,]
Chr4Dxy <- Dxydata[Dxydata$scaffold=="scaffold_4_RagTag",]

Chr4_join <- left_join(Chr4Fst,Chr4Dxy)
Chr4pi <- theta_data[theta_data$Chr=="scaffold_4_RagTag"&theta_data$nSites>16000,]

intervals_df_subset <- intervals_df[intervals_df$Chr=="Chr_4",]
y_top <- max(Chr4Fst$Fst, na.rm = TRUE) * 1.01

y_top <- 0.77

Fst_range <- c(0,0.8)
Dxy_range <- c(0,0.05)
Pi_range <- c(0,0.04)
D_range <- c(-3,1)

a<- ggplot(Chr4_join,aes(start/1000000,Fst))+
  geom_line(color="darkgrey",linewidth=1.5)+
  geom_area(fill='lightgrey')+
  geom_segment(data=intervals_df_subset,aes(x=start/1e6,xend=end/1e6),y=y_top,yend=y_top,color='blue',linewidth=2,inherit.aes=FALSE)+
  labs(y = expression(italic(F)[ST]))+
  xlab("")+
  theme_bw()+
  theme(legend.title=element_blank(),
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),)+
  scale_y_continuous(limits=Fst_range,breaks = Fst_range)
b<- ggplot(Chr4_join,aes(start/1000000,dxy))+
  geom_line(color="darkgrey",linewidth=1.5)+
  geom_area(fill='lightgrey')+
  labs(y = expression(italic(D)[xy]))+
  xlab("")+
  theme_bw()+
  theme(legend.title=element_blank(),
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),)+
  scale_y_continuous(limits=Dxy_range,breaks = Dxy_range)
c<- ggplot(Chr4pi,aes(start/1000000,Mean_pi,group=population,color=population))+
  geom_line(alpha=0.6)+
  labs(y = expression(italic(theta[pi])))+
  xlab("")+
  # xlab("Chromosome coordinates (Mb)")+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position="none",
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),
        legend.text=element_text(size=13),
        plot.margin = unit(c(0, 0, 0, 0), "pt"))+
  scale_color_manual(values=c("red","black"))+
  scale_y_continuous(limits=Pi_range,breaks = Pi_range)
d<- ggplot(Chr4pi,aes(start/1000000,Tajima,group=population,color=population))+
  geom_line(alpha=0.6)+
  labs(y = expression(paste("Tajima's ", italic(D))))+
  xlab("Coordinates (Mb)")+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position="none",
        axis.title=element_text(size=15),
        axis.text=element_text(size=14),
        legend.text=element_text(size=13))+
  scale_color_manual(values=c("red","black"))+
  ylim(-2.7,1)
#scale_y_continuous(limits=D_range,breaks = D_range)




p <- a/b/c/d + 
  plot_layout(ncol = 1) & 
  theme(
    plot.margin = unit(c(5, 10, 0, 10), "pt"),
    legend.position = "none"
  )
p

#ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Chr4_summary_statistics_fst_dxy_pi.png", 
#       plot = p, dpi = 350, width = 10, height = 8, units = "in")
#ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Chr4_summary_statistics_fst_dxy_pi_Tajima.png", 
#       plot = p, dpi = 350, width = 4, height = 5, units = "in")
ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Chr4_summary_statistics_fst_dxy_pi_Tajima_wide.png", 
       plot = p, dpi = 350, width = 8, height = 3.5, units = "in")

##################################################################
#table(!is.element(Chr6Dxy$start,Chr6Fst$start))
# Chr 6 is duplicated in the Dxy dataset for some reason! Probably during the data filtering process
Chr6Fst <- Fstdata[Fstdata$V2=="scaffold_6_RagTag"&Fstdata$V4>16000,]
Chr6Dxy <- Dxydata[Dxydata$scaffold=="scaffold_6_RagTag",]
Chr6Dxy_sort <- Chr6Dxy[order(Chr6Dxy$start),]
Chr6Dxy_sort_nondup <- Chr6Dxy_sort[!duplicated(Chr6Dxy_sort$start),]

Chr6_join <- left_join(Chr6Fst,Chr6Dxy_sort_nondup)

Chr6pi <- theta_data[theta_data$Chr=="scaffold_6_RagTag" &theta_data$nSites>16000,]

intervals_df_subset <- intervals_df[intervals_df$Chr=="Chr_6",]
y_top <- max(Chr6Fst$Fst, na.rm = TRUE) * 1.01

#y_top <- 0.77

Fst_range <- c(0,0.8)
Dxy_range <- c(0,0.05)
Pi_range <- c(0,0.04)
D_range <- c(-3,1)

a<- ggplot(Chr6_join,aes(start/1000000,Fst))+
  geom_line(color="darkgrey",linewidth=1.5)+
  geom_area(fill='lightgrey')+
  geom_segment(data=intervals_df_subset,aes(x=start/1e6,xend=end/1e6),y=y_top,yend=y_top,color='blue',linewidth=2,inherit.aes=FALSE)+
  labs(y = expression(italic(F)[ST]))+
  xlab("")+
  theme_bw()+
  theme(legend.title=element_blank(),
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),)+
scale_y_continuous(limits=Fst_range,breaks = Fst_range)
b<- ggplot(Chr6_join,aes(start/1000000,dxy))+
  geom_line(color="darkgrey",linewidth=1.5)+
  geom_area(fill='lightgrey')+
  labs(y = expression(italic(D)[xy]))+
  xlab("")+
  theme_bw()+
  theme(legend.title=element_blank(),
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),)+
scale_y_continuous(limits=Dxy_range,breaks = Dxy_range)
c<- ggplot(Chr6pi,aes(start/1000000,Mean_pi,group=population,color=population))+
  geom_line(alpha=0.6)+
  labs(y = expression(italic(theta[pi])))+
  xlab("")+
  #xlab("Chromosome coordinates (Mb)")+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position="none",
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),
        legend.text=element_text(size=13),
        plot.margin = unit(c(0, 0, 0, 0), "pt"))+
  scale_color_manual(values=c("red","black"))+
scale_y_continuous(limits=Pi_range,breaks = Pi_range)
d<- ggplot(Chr6pi,aes(start/1000000,Tajima,group=population,color=population))+
  geom_line(alpha=0.6)+
  labs(y = expression(paste("Tajima's ", italic(D))))+
  xlab("Coordinates (Mb)")+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position="none",
        axis.title=element_text(size=15),
        axis.text=element_text(size=14),
        legend.text=element_text(size=13))+
  scale_color_manual(values=c("red","black"))+
  ylim(-2.7,1)
scale_y_continuous(limits=D_range,breaks = D_range)

#p <- a/b/c

#p
#ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Chr6_summary_statistics_fst_dxy_pi.png", 
#       plot = p, dpi = 350, width = 10, height = 8, units = "in")

p <- a/b/c/d + 
  plot_layout(ncol = 1) & 
  theme(
    plot.margin = unit(c(5, 10, 0, 10), "pt"),
    legend.position = "none"
  )
p

#ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Chr6_summary_statistics_fst_dxy_pi_Tajima.png", 
#       plot = p, dpi = 350, width = 4, height = 5, units = "in")

ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Chr6_summary_statistics_fst_dxy_pi_Tajima_wide.png", 
       plot = p, dpi = 350, width = 8, height = 3.5, units = "in")

##################################################################
Chr8Fst <- Fstdata[Fstdata$V2=="scaffold_8_RagTag"& Fstdata$V4>16000,]
Chr8Dxy <- Dxydata[Dxydata$scaffold=="scaffold_8_RagTag",]
Chr8_join <- left_join(Chr8Fst,Chr8Dxy)
Chr8pi <- theta_data[theta_data$Chr=="scaffold_8_RagTag" & theta_data$nSites>16000,]

intervals_df_subset <- intervals_df[intervals_df$Chr=="Chr_8",]
y_top <- max(Chr8Fst$Fst, na.rm = TRUE) * 1.01

y_top <- 0.77

Fst_range <- c(0,0.8)
Dxy_range <- c(0,0.05)
Pi_range <- c(0,0.04)
D_range <- c(-3,1)

a<- ggplot(Chr8_join,aes(start/1000000,Fst))+
  geom_line(color="darkgrey",linewidth=1.5)+
  geom_area(fill='lightgrey')+
  geom_segment(data=intervals_df_subset,aes(x=start/1e6,xend=end/1e6),y=y_top,yend=y_top,color='blue',linewidth=2,inherit.aes=FALSE)+
  labs(y = expression(italic(F)[ST]))+
  xlab("")+
  theme_bw()+
  theme(legend.title=element_blank(),
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),)+
  scale_y_continuous(limits=Fst_range,breaks = Fst_range)
b<- ggplot(Chr8_join,aes(start/1000000,dxy))+
  geom_line(color="darkgrey",linewidth=1.5)+
  geom_area(fill='lightgrey')+
  labs(y = expression(italic(D)[xy]))+
  xlab("")+
  theme_bw()+
  theme(legend.title=element_blank(),
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),)+
  scale_y_continuous(limits=Dxy_range,breaks = Dxy_range)
c<- ggplot(Chr8pi,aes(start/1000000,Mean_pi,group=population,color=population))+
  geom_line(alpha=0.6)+
  labs(y = expression(italic(theta[pi])))+
  xlab("")+
  # xlab("Chromosome coordinates (Mb)")+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position="none",
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),
        legend.text=element_text(size=13),
        plot.margin = unit(c(0, 0, 0, 0), "pt"))+
  scale_color_manual(values=c("red","black"))+
  scale_y_continuous(limits=Pi_range,breaks = Pi_range)
d<- ggplot(Chr8pi,aes(start/1000000,Tajima,group=population,color=population))+
  geom_line(alpha=0.6)+
  labs(y = expression(paste("Tajima's ", italic(D))))+
  xlab("Coordinates (Mb)")+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position="none",
        axis.title=element_text(size=15),
        axis.text=element_text(size=14),
        legend.text=element_text(size=13))+
  scale_color_manual(values=c("red","black"))+
  ylim(-2.7,1)
#scale_y_continuous(limits=D_range,breaks = D_range)




p <- a/b/c/d + 
  plot_layout(ncol = 1) & 
  theme(
    plot.margin = unit(c(5, 10, 0, 10), "pt"),
    legend.position = "none"
  )
p

#ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Chr8_summary_statistics_fst_dxy_pi.png", 
#       plot = p, dpi = 350, width = 10, height = 8, units = "in")
#ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Chr8_summary_statistics_fst_dxy_pi_Tajima.png", 
#       plot = p, dpi = 350, width = 4, height = 5, units = "in")

ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Chr8_summary_statistics_fst_dxy_pi_Tajima_wide.png", 
       plot = p, dpi = 350, width = 8, height = 3.5, units = "in")


##################################################################
Chr9Fst <- Fstdata[Fstdata$V2=="scaffold_9_RagTag" & Fstdata$V4>16000,]
Chr9Dxy <- Dxydata[Dxydata$scaffold=="scaffold_9_RagTag",]
Chr9_join <- left_join(Chr9Fst,Chr9Dxy)
Chr9pi <- theta_data[theta_data$Chr=="scaffold_9_RagTag"& theta_data$nSites>16000,]

intervals_df_subset <- intervals_df[intervals_df$Chr=="Chr_9",]
y_top <- max(Chr15Fst$Fst, na.rm = TRUE) * 1.01

y_top <- 0.77

Fst_range <- c(0,0.8)
Dxy_range <- c(0,0.05)
Pi_range <- c(0,0.04)
D_range <- c(-3,1)

a<- ggplot(Chr9_join,aes(start/1000000,Fst))+
  geom_line(color="darkgrey",linewidth=1.5)+
  geom_area(fill='lightgrey')+
  geom_segment(data=intervals_df_subset,aes(x=start/1e6,xend=end/1e6),y=y_top,yend=y_top,color='blue',linewidth=2,inherit.aes=FALSE)+
  labs(y = expression(italic(F)[ST]))+
  xlab("")+
  theme_bw()+
  theme(legend.title=element_blank(),
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),)+
  scale_y_continuous(limits=Fst_range,breaks = Fst_range)
b<- ggplot(Chr9_join,aes(start/1000000,dxy))+
  geom_line(color="darkgrey",linewidth=1.5)+
  geom_area(fill='lightgrey')+
  labs(y = expression(italic(D)[xy]))+
  xlab("")+
  theme_bw()+
  theme(legend.title=element_blank(),
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),)+
  scale_y_continuous(limits=Dxy_range,breaks = Dxy_range)
c<- ggplot(Chr9pi,aes(start/1000000,Mean_pi,group=population,color=population))+
  geom_line(alpha=0.6)+
  labs(y = expression(italic(theta[pi])))+
  xlab("")+
  # xlab("Chromosome coordinates (Mb)")+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position="none",
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),
        legend.text=element_text(size=13),
        plot.margin = unit(c(0, 0, 0, 0), "pt"))+
  scale_color_manual(values=c("red","black"))+
  scale_y_continuous(limits=Pi_range,breaks = Pi_range)
d<- ggplot(Chr9pi,aes(start/1000000,Tajima,group=population,color=population))+
  geom_line(alpha=0.6)+
  labs(y = expression(paste("Tajima's ", italic(D))))+
  xlab("Coordinates (Mb)")+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position="none",
        axis.title=element_text(size=15),
        axis.text=element_text(size=14),
        legend.text=element_text(size=13))+
  scale_color_manual(values=c("red","black"))+
  ylim(-2.7,1)
#scale_y_continuous(limits=D_range,breaks = D_range)




p <- a/b/c/d + 
  plot_layout(ncol = 1) & 
  theme(
    plot.margin = unit(c(5, 10, 0, 10), "pt"),
    legend.position = "none"
  )
p

#ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Chr9_summary_statistics_fst_dxy_pi.png", 
#       plot = p, dpi = 350, width = 10, height = 8, units = "in")
ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Chr9_summary_statistics_fst_dxy_pi_Tajima.png", 
       plot = p, dpi = 350, width = 4, height = 5, units = "in")





##################################################################
Chr10Fst <- Fstdata[Fstdata$V2=="scaffold_10_RagTag"&Fstdata$V4>16000,]
Chr10Dxy <- Dxydata[Dxydata$scaffold=="scaffold_10_RagTag",]
Chr10_join <- left_join(Chr10Fst,Chr10Dxy)
Chr10pi <- theta_data[theta_data$Chr=="scaffold_10_RagTag"&theta_data$nSites>16000,]

intervals_df_subset <- intervals_df[intervals_df$Chr=="Chr_10",]

y_top <- 0.77

Fst_range <- c(0,0.8)
Dxy_range <- c(0,0.05)
Pi_range <- c(0,0.04)
D_range <- c(-3,1)

a<- ggplot(Chr10_join,aes(start/1000000,Fst))+
  geom_line(color="darkgrey",linewidth=1.5)+
  geom_area(fill='lightgrey')+
  geom_segment(data=intervals_df_subset,aes(x=start/1e6,xend=end/1e6),y=y_top,yend=y_top,color='blue',linewidth=2,inherit.aes=FALSE)+
  labs(y = expression(italic(F)[ST]))+
  xlab("")+
  theme_bw()+
  theme(legend.title=element_blank(),
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),)+
  scale_y_continuous(limits=Fst_range,breaks = Fst_range)
b<- ggplot(Chr10_join,aes(start/1000000,dxy))+
  geom_line(color="darkgrey",linewidth=1.5)+
  geom_area(fill='lightgrey')+
  labs(y = expression(italic(D)[xy]))+
  xlab("")+
  theme_bw()+
  theme(legend.title=element_blank(),
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),)+
  scale_y_continuous(limits=Dxy_range,breaks = Dxy_range)
c<- ggplot(Chr10pi,aes(start/1000000,Mean_pi,group=population,color=population))+
  geom_line(alpha=0.6)+
  labs(y = expression(italic(theta[pi])))+
  xlab("")+
  # xlab("Chromosome coordinates (Mb)")+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position="none",
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),
        legend.text=element_text(size=13),
        plot.margin = unit(c(0, 0, 0, 0), "pt"))+
  scale_color_manual(values=c("red","black"))+
  scale_y_continuous(limits=Pi_range,breaks = Pi_range)
d<- ggplot(Chr10pi,aes(start/1000000,Tajima,group=population,color=population))+
  geom_line(alpha=0.6)+
  labs(y = expression(paste("Tajima's ", italic(D))))+
  xlab("Coordinates (Mb)")+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position="none",
        axis.title=element_text(size=15),
        axis.text=element_text(size=14),
        legend.text=element_text(size=13))+
  scale_color_manual(values=c("red","black"))+
  ylim(-2.7,1)


p <- a/b/c/d + 
  plot_layout(ncol = 1) & 
  theme(
    plot.margin = unit(c(5, 10, 0, 10), "pt"),
    legend.position = "none"
  )
p

#ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Chr10_summary_statistics_fst_dxy_pi.png", 
#       plot = p, dpi = 350, width = 10, height = 8, units = "in")
ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Chr10_summary_statistics_fst_dxy_pi_Tajima.png", 
       plot = p, dpi = 350, width = 4, height = 5, units = "in")




##################################################################
Chr13Fst <- Fstdata[Fstdata$V2=="scaffold_13_RagTag"&Fstdata$V4>16000,]
Chr13Dxy <- Dxydata[Dxydata$scaffold=="scaffold_13_RagTag",]
Chr13_join <- left_join(Chr13Fst,Chr13Dxy)
Chr13pi <- theta_data[theta_data$Chr=="scaffold_13_RagTag"&theta_data$nSites>16000,]

intervals_df_subset <- intervals_df[intervals_df$Chr=="Chr_13",]
y_top <- max(Chr15Fst$Fst, na.rm = TRUE) * 1.01

y_top <- 0.77

Fst_range <- c(0,0.8)
Dxy_range <- c(0,0.05)
Pi_range <- c(0,0.04)
D_range <- c(-3,1)

a<- ggplot(Chr13_join,aes(start/1000000,Fst))+
  geom_line(color="darkgrey",linewidth=1.5)+
  geom_area(fill='lightgrey')+
  geom_segment(data=intervals_df_subset,aes(x=start/1e6,xend=end/1e6),y=y_top,yend=y_top,color='blue',linewidth=2,inherit.aes=FALSE)+
  labs(y = expression(italic(F)[ST]))+
  xlab("")+
  theme_bw()+
  theme(legend.title=element_blank(),
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),)+
  scale_y_continuous(limits=Fst_range,breaks = Fst_range)
b<- ggplot(Chr13_join,aes(start/1000000,dxy))+
  geom_line(color="darkgrey",linewidth=1.5)+
  geom_area(fill='lightgrey')+
  labs(y = expression(italic(D)[xy]))+
  xlab("")+
  theme_bw()+
  theme(legend.title=element_blank(),
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),)+
  scale_y_continuous(limits=Dxy_range,breaks = Dxy_range)
c<- ggplot(Chr13pi,aes(start/1000000,Mean_pi,group=population,color=population))+
  geom_line(alpha=0.6)+
  labs(y = expression(italic(theta[pi])))+
  xlab("")+
  # xlab("Chromosome coordinates (Mb)")+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position="none",
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),
        legend.text=element_text(size=13),
        plot.margin = unit(c(0, 0, 0, 0), "pt"))+
  scale_color_manual(values=c("red","black"))+
  scale_y_continuous(limits=Pi_range,breaks = Pi_range)
d<- ggplot(Chr13pi,aes(start/1000000,Tajima,group=population,color=population))+
  geom_line(alpha=0.6)+
  labs(y = expression(paste("Tajima's ", italic(D))))+
  xlab("Coordinates (Mb)")+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position="none",
        axis.title=element_text(size=15),
        axis.text=element_text(size=14),
        legend.text=element_text(size=13))+
  scale_color_manual(values=c("red","black"))+
  ylim(-2.7,1)
#scale_y_continuous(limits=D_range,breaks = D_range)




p <- a/b/c/d + 
  plot_layout(ncol = 1) & 
  theme(
    plot.margin = unit(c(5, 10, 0, 10), "pt"),
    legend.position = "none"
  )
p

#ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Chr13_summary_statistics_fst_dxy_pi.png", 
#       plot = p, dpi = 350, width = 10, height = 8, units = "in")
ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Chr13_summary_statistics_fst_dxy_pi_Tajima.png", 
       plot = p, dpi = 350, width = 4, height = 5, units = "in")

##################################################################
Chr14Fst <- Fstdata[Fstdata$V2=="scaffold_14_RagTag"&Fstdata$V4>16000,]
Chr14Dxy <- Dxydata[Dxydata$scaffold=="scaffold_14_RagTag",]
Chr14_join <- left_join(Chr14Fst,Chr14Dxy)
Chr14pi <- theta_data[theta_data$Chr=="scaffold_14_RagTag"&theta_data$nSites>16000,]

intervals_df_subset <- intervals_df[intervals_df$Chr=="Chr_14",]
y_top <- max(Chr14Fst$Fst, na.rm = TRUE) * 1.01

y_top <- 0.77

Fst_range <- c(0,0.8)
Dxy_range <- c(0,0.05)
Pi_range <- c(0,0.04)
D_range <- c(-3,1)

a<- ggplot(Chr14_join,aes(start/1000000,Fst))+
  geom_line(color="darkgrey",linewidth=1.5)+
  geom_area(fill='lightgrey')+
  geom_segment(data=intervals_df_subset,aes(x=start/1e6,xend=end/1e6),y=y_top,yend=y_top,color='blue',linewidth=2,inherit.aes=FALSE)+
  labs(y = expression(italic(F)[ST]))+
  xlab("")+
  theme_bw()+
  theme(legend.title=element_blank(),
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),)+
 scale_y_continuous(limits=Fst_range,breaks = Fst_range)
b<- ggplot(Chr14_join,aes(start/1000000,dxy))+
  geom_line(color="darkgrey",linewidth=1.5)+
  geom_area(fill='lightgrey')+
  labs(y = expression(italic(D)[xy]))+
  xlab("")+
  theme_bw()+
  theme(legend.title=element_blank(),
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),)+
scale_y_continuous(limits=Dxy_range,breaks = Dxy_range)
c<- ggplot(Chr14pi,aes(start/1000000,Mean_pi,group=population,color=population))+
  geom_line(alpha=0.6)+
  labs(y = expression(italic(theta[pi])))+
  xlab("")+
 # xlab("Chromosome coordinates (Mb)")+
  theme_bw()+
  theme(legend.title=element_blank(),
       # legend.position="none",
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),
        legend.text=element_text(size=13))+
  scale_color_manual(values=c("red","black"))+
scale_y_continuous(limits=Pi_range,breaks = Pi_range)
d<- ggplot(Chr14pi,aes(start/1000000,Tajima,group=population,color=population))+
  geom_line(alpha=0.6)+
  labs(y = expression(paste("Tajima's ", italic(D))))+
  xlab("Coordinates (Mb)")+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position="none",
        axis.title=element_text(size=15),
        axis.text=element_text(size=14),
        legend.text=element_text(size=13))+
  scale_color_manual(values=c("red","black"))+
  ylim(-2.7,1)
#scale_y_continuous(limits=D_range,breaks = D_range)


#p <- a/b/c 

#p


#ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Chr14_summary_statistics_fst_dxy_pi.png", 
#       plot = p, dpi = 350, width = 10, height = 8, units = "in")


p <- a/b/c/d + 
  plot_layout(ncol = 1) & 
  theme(
    plot.margin = unit(c(5, 10, 0, 10), "pt"),
    legend.position = "none"
  )
p

ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Chr14_summary_statistics_fst_dxy_pi_Tajima.png", 
       plot = p, dpi = 350, width = 4, height = 5, units = "in")


#ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Chr14_summary_statistics_fst_dxy_pi_Tajima_wide.png", 
#       plot = p, dpi = 350, width = 7, height = 7, units = "in")


##################################################################
Chr15Fst <- Fstdata[Fstdata$V2=="scaffold_15_RagTag"&Fstdata$V4>16000,]
Chr15Dxy <- Dxydata[Dxydata$scaffold=="scaffold_15_RagTag",]
Chr15_join <- left_join(Chr15Fst,Chr15Dxy)
Chr15pi <- theta_data[theta_data$Chr=="scaffold_15_RagTag"&theta_data$nSites>16000,]

intervals_df_subset <- intervals_df[intervals_df$Chr=="Chr_15",]
y_top <- max(Chr15Fst$Fst, na.rm = TRUE) * 1.01

y_top <- 0.77

Fst_range <- c(0,0.8)
Dxy_range <- c(0,0.05)
Pi_range <- c(0,0.04)
D_range <- c(-3,1)

a<- ggplot(Chr15_join,aes(start/1000000,Fst))+
  geom_line(color="darkgrey",linewidth=1.5)+
  geom_area(fill='lightgrey')+
  geom_segment(data=intervals_df_subset,aes(x=start/1e6,xend=end/1e6),y=y_top,yend=y_top,color='blue',linewidth=2,inherit.aes=FALSE)+
  labs(y = expression(italic(F)[ST]))+
  xlab("")+
  theme_bw()+
  theme(legend.title=element_blank(),
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),)+
  scale_y_continuous(limits=Fst_range,breaks = Fst_range)
b<- ggplot(Chr15_join,aes(start/1000000,dxy))+
  geom_line(color="darkgrey",linewidth=1.5)+
  geom_area(fill='lightgrey')+
  labs(y = expression(italic(D)[xy]))+
  xlab("")+
  theme_bw()+
  theme(legend.title=element_blank(),
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),)+
  scale_y_continuous(limits=Dxy_range,breaks = Dxy_range)
c<- ggplot(Chr15pi,aes(start/1000000,Mean_pi,group=population,color=population))+
  geom_line(alpha=0.6)+
  labs(y = expression(italic(theta[pi])))+
  xlab("")+
  # xlab("Chromosome coordinates (Mb)")+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position="none",
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),
        legend.text=element_text(size=13),
        plot.margin = unit(c(0, 0, 0, 0), "pt"))+
  scale_color_manual(values=c("red","black"))+
  scale_y_continuous(limits=Pi_range,breaks = Pi_range)
d<- ggplot(Chr15pi,aes(start/1000000,Tajima,group=population,color=population))+
  geom_line(alpha=0.6)+
  labs(y = expression(paste("Tajima's ", italic(D))))+
  xlab("Coordinates (Mb)")+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position="none",
        axis.title=element_text(size=15),
        axis.text=element_text(size=14),
        legend.text=element_text(size=13))+
  scale_color_manual(values=c("red","black"))+
  ylim(-2.7,1)
#scale_y_continuous(limits=D_range,breaks = D_range)




p <- a/b/c/d + 
  plot_layout(ncol = 1) & 
  theme(
    plot.margin = unit(c(5, 10, 0, 10), "pt"),
    legend.position = "none"
  )
p

#ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Chr15_summary_statistics_fst_dxy_pi.png", 
#       plot = p, dpi = 350, width = 10, height = 8, units = "in")
#ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Chr15_summary_statistics_fst_dxy_pi_Tajima.png", 
#       plot = p, dpi = 350, width = 4, height = 5, units = "in")



##################################################################
Chr17Fst <- Fstdata[Fstdata$V2=="scaffold_17_RagTag"&Fstdata$V4>16000,]
Chr17Dxy <- Dxydata[Dxydata$scaffold=="scaffold_17_RagTag",]
Chr17_join <- left_join(Chr17Fst,Chr17Dxy)
Chr17pi <- theta_data[theta_data$Chr=="scaffold_17_RagTag"&theta_data$nSites>16000,]


intervals_df_subset <- intervals_df[intervals_df$Chr=="Chr_17",]
y_top <- max(Chr17Fst$Fst, na.rm = TRUE) * 1.01

y_top <- 0.77

Fst_range <- c(0,0.8)
Dxy_range <- c(0,0.05)
Pi_range <- c(0,0.04)
D_range <- c(-3,1)

a<- ggplot(Chr17_join,aes(start/1000000,Fst))+
  geom_line(color="darkgrey",linewidth=1.5)+
  geom_area(fill='lightgrey')+
  geom_segment(data=intervals_df_subset,aes(x=start/1e6,xend=end/1e6),y=y_top,yend=y_top,color='blue',linewidth=2,inherit.aes=FALSE)+
  labs(y = expression(italic(F)[ST]))+
  xlab("")+
  theme_bw()+
  theme(legend.title=element_blank(),
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),)+
  scale_y_continuous(limits=Fst_range,breaks = Fst_range)
b<- ggplot(Chr17_join,aes(start/1000000,dxy))+
  geom_line(color="darkgrey",linewidth=1.5)+
  geom_area(fill='lightgrey')+
  labs(y = expression(italic(D)[xy]))+
  xlab("")+
  theme_bw()+
  theme(legend.title=element_blank(),
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),)+
  scale_y_continuous(limits=Dxy_range,breaks = Dxy_range)
c<- ggplot(Chr17pi,aes(start/1000000,Mean_pi,group=population,color=population))+
  geom_line(alpha=0.6)+
  labs(y = expression(italic(theta[pi])))+
  xlab("")+
  # xlab("Chromosome coordinates (Mb)")+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position="none",
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),
        legend.text=element_text(size=13),
        plot.margin = unit(c(0, 0, 0, 0), "pt"))+
  scale_color_manual(values=c("red","black"))+
  scale_y_continuous(limits=Pi_range,breaks = Pi_range)
d<- ggplot(Chr17pi,aes(start/1000000,Tajima,group=population,color=population))+
  geom_line(alpha=0.6)+
  labs(y = expression(paste("Tajima's ", italic(D))))+
  xlab("Coordinates (Mb)")+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position="none",
        axis.title=element_text(size=15),
        axis.text=element_text(size=14),
        legend.text=element_text(size=13))+
  scale_color_manual(values=c("red","black"))+
  ylim(-2.7,1)
#scale_y_continuous(limits=D_range,breaks = D_range)




p <- a/b/c/d + 
  plot_layout(ncol = 1) & 
  theme(
    plot.margin = unit(c(5, 10, 0, 10), "pt"),
    legend.position = "none"
  )
p

#ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Chr17_summary_statistics_fst_dxy_pi.png", 
#       plot = p, dpi = 350, width = 10, height = 8, units = "in")
ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Chr17_summary_statistics_fst_dxy_pi_Tajima.png", 
       plot = p, dpi = 350, width = 4, height = 5, units = "in")
#ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Chr17_summary_statistics_fst_dxy_pi_Tajima.png", 
#       plot = p, dpi = 350, width = 10, height = 8, units = "in")


##################################################################
Chr18Fst <- Fstdata[Fstdata$V2=="scaffold_18_RagTag"&Fstdata$V4>16000,]
Chr18Dxy <- Dxydata[Dxydata$scaffold=="scaffold_18_RagTag",]
Chr18_join <- left_join(Chr18Fst,Chr18Dxy)
Chr18pi <- theta_data[theta_data$Chr=="scaffold_18_RagTag"&theta_data$nSites>16000,]

intervals_df_subset <- intervals_df[intervals_df$Chr=="Chr_18",]
y_top <- max(Chr18Fst$Fst, na.rm = TRUE) * 1.01

y_top <- 0.77

Fst_range <- c(0,0.8)
Dxy_range <- c(0,0.05)
Pi_range <- c(0,0.04)
D_range <- c(-3,1)

a<- ggplot(Chr18_join,aes(start/1000000,Fst))+
  geom_line(color="darkgrey",linewidth=1.5)+
  geom_area(fill='lightgrey')+
  geom_segment(data=intervals_df_subset,aes(x=start/1e6,xend=end/1e6),y=y_top,yend=y_top,color='blue',linewidth=2,inherit.aes=FALSE)+
  labs(y = expression(italic(F)[ST]))+
  xlab("")+
  theme_bw()+
  theme(legend.title=element_blank(),
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),)+
  scale_y_continuous(limits=Fst_range,breaks = Fst_range)
b<- ggplot(Chr18_join,aes(start/1000000,dxy))+
  geom_line(color="darkgrey",linewidth=1.5)+
  geom_area(fill='lightgrey')+
  labs(y = expression(italic(D)[xy]))+
  xlab("")+
  theme_bw()+
  theme(legend.title=element_blank(),
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),)+
  scale_y_continuous(limits=Dxy_range,breaks = Dxy_range)
c<- ggplot(Chr18pi,aes(start/1000000,Mean_pi,group=population,color=population))+
  geom_line(alpha=0.6)+
  labs(y = expression(italic(theta[pi])))+
  xlab("")+
  # xlab("Chromosome coordinates (Mb)")+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position="none",
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),
        legend.text=element_text(size=13),
        plot.margin = unit(c(0, 0, 0, 0), "pt"))+
  scale_color_manual(values=c("red","black"))+
  scale_y_continuous(limits=Pi_range,breaks = Pi_range)
d<- ggplot(Chr18pi,aes(start/1000000,Tajima,group=population,color=population))+
  geom_line(alpha=0.6)+
  labs(y = expression(paste("Tajima's ", italic(D))))+
  xlab("Coordinates (Mb)")+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position="none",
        axis.title=element_text(size=15),
        axis.text=element_text(size=14),
        legend.text=element_text(size=13))+
  scale_color_manual(values=c("red","black"))+
  ylim(-2.7,1)
#scale_y_continuous(limits=D_range,breaks = D_range)




p <- a/b/c/d + 
  plot_layout(ncol = 1) & 
  theme(
    plot.margin = unit(c(5, 10, 0, 10), "pt"),
    legend.position = "none"
  )
p

#ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Chr18_summary_statistics_fst_dxy_pi.png", 
#       plot = p, dpi = 350, width = 10, height = 8, units = "in")
ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Chr18_summary_statistics_fst_dxy_pi_Tajima.png", 
       plot = p, dpi = 350, width = 4, height = 5, units = "in")

##################################################################
pri15Fst <- Fstdata[Fstdata$V2=="pri#ptg000015l" & Fstdata$V4>16000,]
pri15Dxy <- Dxydata[Dxydata$scaffold=="pri#ptg000015l",]
pri15_join <- left_join(pri15Fst,pri15Dxy)
pri15pi <- theta_data[theta_data$Chr=="pri#ptg000015l"&theta_data$nSites>16000,]
intervals_df_subset <- intervals_df[intervals_df$Chr=="pri#ptg000015l",]
y_top <- max(pri15Fst$Fst, na.rm = TRUE) * 1.01


y_top <- 0.77

Fst_range <- c(0,0.8)
Dxy_range <- c(0,0.04)
Pi_range <- c(0,0.04)
D_range <- c(-3,1)

a<- ggplot(pri15_join,aes(start/1000000,Fst))+
  geom_line(color="darkgrey",linewidth=1.5)+
  geom_area(fill='lightgrey')+
  geom_segment(data=intervals_df_subset,aes(x=start/1e6,xend=end/1e6),y=y_top,yend=y_top,color='blue',linewidth=2,inherit.aes=FALSE)+
  labs(y = expression(italic(F)[ST]))+
  xlab("")+
  theme_bw()+
  theme(legend.title=element_blank(),
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),)+
  scale_y_continuous(limits=Fst_range,breaks = Fst_range)
b<- ggplot(pri15_join,aes(start/1000000,dxy))+
  geom_line(color="darkgrey",linewidth=1.5)+
  geom_area(fill='lightgrey')+
  labs(y = expression(italic(D)[xy]))+
  xlab("")+
  theme_bw()+
  theme(legend.title=element_blank(),
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),)+
  scale_y_continuous(limits=Dxy_range,breaks = Dxy_range)
c<- ggplot(pri15pi,aes(start/1000000,Mean_pi,group=population,color=population))+
  geom_line(alpha=0.6)+
  labs(y = expression(italic(theta[pi])))+
  xlab("")+
  # xlab("Chromosome coordinates (Mb)")+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position="none",
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),
        legend.text=element_text(size=13),
        plot.margin = unit(c(0, 0, 0, 0), "pt"))+
  scale_color_manual(values=c("red","black"))+
  scale_y_continuous(limits=Pi_range,breaks = Pi_range)
d<- ggplot(pri15pi,aes(start/1000000,Tajima,group=population,color=population))+
  geom_line(alpha=0.6)+
  labs(y = expression(paste("Tajima's ", italic(D))))+
  xlab("Coordinates (Mb)")+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position="none",
        axis.title=element_text(size=15),
        axis.text=element_text(size=14),
        legend.text=element_text(size=13))+
  scale_color_manual(values=c("red","black"))+
  ylim(-2.7,1)
#scale_y_continuous(limits=D_range,breaks = D_range)




p <- a/b/c/d + 
  plot_layout(ncol = 1) & 
  theme(
    plot.margin = unit(c(5, 10, 0, 10), "pt"),
    legend.position = "none"
  )
p

#ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/pri15_summary_statistics_fst_dxy_pi.png", 
#       plot = p, dpi = 350, width = 10, height = 8, units = "in")
#ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/pri15_summary_statistics_fst_dxy_pi_Tajima.png", 
#       plot = p, dpi = 350, width = 4, height = 5, units = "in")

ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/pri15_summary_statistics_fst_dxy_pi_Tajima_wide.png", 
       plot = p, dpi = 350, width = 7, height = 7, units = "in")

##################################################################
pri45Fst <- Fstdata[Fstdata$V2=="pri#ptg000045l" & Fstdata$V4>16000,]
pri45Dxy <- Dxydata[Dxydata$scaffold=="pri#ptg000045l",]

pri45_join <- left_join(pri45Fst,pri45Dxy)
pri45pi <- theta_data[theta_data$Chr=="pri#ptg000045l"&theta_data$nSites>16000,]

intervals_df_subset <- intervals_df[intervals_df$Chr=="pri#ptg000045l",]
y_top <- max(pri45Fst$Fst, na.rm = TRUE) * 1.01

y_top <- 0.77

Fst_range <- c(0,0.8)
Dxy_range <- c(0,0.05)
Pi_range <- c(0,0.04)
D_range <- c(-3,1)

a<- ggplot(pri45_join,aes(start/1000000,Fst))+
  geom_line(color="darkgrey",linewidth=1.5)+
  geom_area(fill='lightgrey')+
  geom_segment(data=intervals_df_subset,aes(x=start/1e6,xend=end/1e6),y=y_top,yend=y_top,color='blue',linewidth=2,inherit.aes=FALSE)+
  labs(y = expression(italic(F)[ST]))+
  xlab("")+
  theme_bw()+
  theme(legend.title=element_blank(),
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),)+
  scale_y_continuous(limits=Fst_range,breaks = Fst_range)
b<- ggplot(pri45_join,aes(start/1000000,dxy))+
  geom_line(color="darkgrey",linewidth=1.5)+
  geom_area(fill='lightgrey')+
  labs(y = expression(italic(D)[xy]))+
  xlab("")+
  theme_bw()+
  theme(legend.title=element_blank(),
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),)+
  scale_y_continuous(limits=Dxy_range,breaks = Dxy_range)
c<- ggplot(pri45pi,aes(start/1000000,Mean_pi,group=population,color=population))+
  geom_line(alpha=0.6)+
  labs(y = expression(italic(theta[pi])))+
  xlab("")+
 # xlab("Chromosome coordinates (Mb)")+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position="none",
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),
        legend.text=element_text(size=13),
        plot.margin = unit(c(0, 0, 0, 0), "pt"))+
  scale_color_manual(values=c("red","black"))+
  scale_y_continuous(limits=Pi_range,breaks = Pi_range)
d<- ggplot(pri45pi,aes(start/1000000,Tajima,group=population,color=population))+
  geom_line(alpha=0.6)+
  labs(y = expression(paste("Tajima's ", italic(D))))+
  xlab("Coordinates (Mb)")+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position="none",
        axis.title=element_text(size=15),
        axis.text=element_text(size=14),
        legend.text=element_text(size=13))+
  scale_color_manual(values=c("red","black"))+
  ylim(-2.7,1)
  #scale_y_continuous(limits=D_range,breaks = D_range)




p <- a/b/c/d + 
  plot_layout(ncol = 1) & 
  theme(
    plot.margin = unit(c(5, 10, 0, 10), "pt"),
    legend.position = "none"
  )
p

#ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/pri45_summary_statistics_fst_dxy_pi.png", 
#       plot = p, dpi = 350, width = 10, height = 8, units = "in")
ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/pri45_summary_statistics_fst_dxy_pi_Tajima.png", 
       plot = p, dpi = 350, width = 4, height = 5, units = "in")


##################################################################
neoZFst <- Fstdata[Fstdata$V2=="scaffold_2_RagTag" & Fstdata$V4>16000,]
neoZDxy <- Dxydata[Dxydata$scaffold=="scaffold_2_RagTag",]

neoZ_join <- left_join(neoZFst,neoZDxy)
neoZpi <- theta_data[theta_data$Chr=="scaffold_2_RagTag"&theta_data$nSites>16000,]
y_top <- max(neoZFst$Fst, na.rm = TRUE) * 1.05
# Plot intervals for neo-Z
intervals_neoZ <- data.frame(Chr=c("Neo_Z","Neo_Z","Neo_Z"),start=c(0,74.4,111.6),end=c(74.4,111.6,128.9))
intervals_neoZ$Chr <- as.factor(intervals_neoZ$Chr)

a<- ggplot(neoZ_join,aes(start/1000000,Fst))+
  geom_line()+
  geom_area(fill='lightgrey')+
  ylim(0,0.24)+
  geom_segment(data = intervals_neoZ,aes(x=start,xend=end),y = y_top, yend = y_top,color = 'black', linewidth = 2.8,inherit.aes=FALSE) +
  geom_segment(data=intervals_neoZ,aes(x=start,xend=end),y=y_top,yend=y_top,color=c('black','white','#857E13'),linewidth=2,inherit.aes=FALSE)+
  labs(y = expression(italic(F)[ST]))+
  xlab("")+
  theme_bw()+
  theme(legend.title=element_blank(),
        axis.title=element_text(size=15),
        axis.text=element_text(size=14),)
b<- ggplot(neoZ_join,aes(start/1000000,dxy))+
  geom_line()+
  geom_area(fill='lightgrey')+
  labs(y = expression(italic(D)[xy]))+
  xlab("")+
  theme_bw()+
  theme(legend.title=element_blank(),
        axis.title=element_text(size=15),
        axis.text=element_text(size=14),)
c<- ggplot(neoZpi,aes(start/1000000,Mean_pi,group=population,color=population))+
  geom_line()+
  labs(y = expression(italic(theta[pi])))+
  xlab("")+
  # xlab("Chromosome coordinates (Mb)")+
  theme_bw()+
  theme(legend.title=element_blank(),
        axis.title=element_text(size=15),
        axis.text=element_text(size=14),
        legend.text=element_text(size=13))+
  scale_color_manual(values=c("red","black"))
d<- ggplot(neoZpi,aes(start/1000000,Tajima,group=population,color=population))+
  geom_line()+
  labs(y = expression(paste("Tajima's ", italic(D))))+
  xlab("neo-Z coordinates (Mb)")+
  theme_bw()+
  theme(legend.title=element_blank(),
        axis.title=element_text(size=15),
        axis.text=element_text(size=14),
        legend.text=element_text(size=13))+
  scale_color_manual(values=c("red","black"))




p <- a/b/c/d
p

ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/neoZ_summary_statistics_fst_dxy_pi_Tajima.png", 
       plot = p, dpi = 350, width = 10, height = 8, units = "in")




# Chromosomes WITHOUT MOR outliers!!!!!

##################################################################
Chr1Fst <- Fstdata[Fstdata$V2=="scaffold_1_RagTag"&Fstdata$V4>16000,]
Chr1Dxy <- Dxydata[Dxydata$scaffold=="scaffold_1_RagTag",]

Chr1_join <- left_join(Chr1Fst,Chr1Dxy)
Chr1pi <- theta_data[theta_data$Chr=="scaffold_1_RagTag"&theta_data$nSites>16000,]


y_top <- max(Chr1Fst$Fst, na.rm = TRUE) * 1.01

#y_top <- 0.77

Fst_range <- c(0,0.2)
Dxy_range <- c(0,0.05)
Pi_range <- c(0,0.04)
D_range <- c(-3,1)

a<- ggplot(Chr1_join,aes(start/1000000,Fst))+
  geom_line(color="darkgrey",linewidth=1.5)+
  geom_area(fill='lightgrey')+
  labs(y = expression(italic(F)[ST]))+
  xlab("")+
  theme_bw()+
  theme(legend.title=element_blank(),
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),)+
  scale_y_continuous(limits=Fst_range,breaks = Fst_range)
b<- ggplot(Chr1_join,aes(start/1000000,dxy))+
  geom_line(color="darkgrey",linewidth=1.5)+
  geom_area(fill='lightgrey')+
  labs(y = expression(italic(D)[xy]))+
  xlab("")+
  theme_bw()+
  theme(legend.title=element_blank(),
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),)+
  scale_y_continuous(limits=Dxy_range,breaks = Dxy_range)
c<- ggplot(Chr1pi,aes(start/1000000,Mean_pi,group=population,color=population))+
  geom_line(alpha=0.6)+
  labs(y = expression(italic(theta[pi])))+
  xlab("")+
  # xlab("Chromosome coordinates (Mb)")+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position="none",
        axis.title=element_text(size=15),
        axis.text.x=element_blank(),
        axis.text=element_text(size=14),
        legend.text=element_text(size=13),
        plot.margin = unit(c(0, 0, 0, 0), "pt"))+
  scale_color_manual(values=c("red","black"))+
  scale_y_continuous(limits=Pi_range,breaks = Pi_range)
d<- ggplot(Chr1pi,aes(start/1000000,Tajima,group=population,color=population))+
  geom_line(alpha=0.6)+
  labs(y = expression(paste("Tajima's ", italic(D))))+
  xlab("Coordinates (Mb)")+
  theme_bw()+
  theme(legend.title=element_blank(),
        legend.position="none",
        axis.title=element_text(size=15),
        axis.text=element_text(size=14),
        legend.text=element_text(size=13))+
  scale_color_manual(values=c("red","black"))+
  ylim(-2.7,1)
#scale_y_continuous(limits=D_range,breaks = D_range)




p <- a/b/c/d + 
  plot_layout(ncol = 1) & 
  theme(
    plot.margin = unit(c(5, 10, 0, 10), "pt"),
    legend.position = "none"
  )
p

#ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Chr1_summary_statistics_fst_dxy_pi.png", 
#       plot = p, dpi = 350, width = 10, height = 8, units = "in")
ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Chr1_summary_statistics_fst_dxy_pi_Tajima.png", 
       plot = p, dpi = 350, width = 4, height = 5, units = "in")

