# Calculate relatedness from all autosomal scaffolds <200kb
# January 2025
# S MacRae Orzechowski

library(ggplot2)
library(stringr)
library(ggrepel)
#install.packages("ggrepel")
setwd("~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/ngsRelate")

data <- read.table('autosomes_ngsRelate_16threads.txt',header=T)
samples<-read.table('allbamspath.txt',sep="/")
samples$samples <- str_remove_all(samples$V12,'_bwa_merge_dedup_clip_sort_fixmate_realigned_final.bam|_bwa_dedup_clip_sort_fixmate_realigned_final.bam')


pairs_to_label_firstorder <- data[data$KING >0.17,]
pairs_to_label_secorder <- data[data$KING < 0.17 & data$KING > 0.08,]
pairs_to_label_thirdorder <- data[data$KING < 0.08 & data$KING > 0.04,]

# Well, these values do seem rather inflated, given the third order frequency, don't you think?
# I'm not sure why -- I shall have to check! Well, it's not necessarily inflated, some individuals
# with little data, like Binya_366910_M are inferred to be third order individuals, but that's just data limitations
# Also, it's not 73 individuals, it's 73 pairs, often involving the same individual, so it's not as bad as I thought!


## Redid this! The first time I was using the wrong IDs because of 0-based vs 1-based numbering in the data file!

pairs_to_label_firstorder$names <- c("Mallee_313_M_314_M","Moon_708_F_024_M","Nom_100_M_101_F","Pill_110_M_115_M","Pill_111_M_112_M",
                                     "Walch_120_F_154_M","Walch_121_F_154_M")

pairs_to_label_secorder$names <- c("Berth_157_M_162_F","Mallee_321_F_324_M","Walch_119_M_Walch_122_M","Walch_120_F_121_F")

pairs_to_label <- rbind(pairs_to_label_firstorder,pairs_to_label_secorder)



## First order
# Individuals 24 and 25, which correspond to 25 and 26 in the sample names because of 0-based vs 1-based
samples$samples[25]
samples$samples[26]

# Individuals 31 and 35, which correspond to 32 and 36 in the sample names
samples$samples[32]
samples$samples[36]

# Individuals 46 and 47, which correspond to 47 and 48 in the sample names
samples$samples[47]
samples$samples[48]

# Individuals 51 and 56, which correspond to 52 and 57 in the sample names
samples$samples[52]
samples$samples[57]

# Individuals 52 and 53, which correspond to 53 and 54 in the sample names
samples$samples[53]
samples$samples[54]

# Individuals 64 and 69, which correspond to 65 and 70 in the sample names
samples$samples[65]
samples$samples[70]

# Individuals 65 and 69, which correspond to 66 and 70 in the sample names
samples$samples[66]
samples$samples[70]

## Second order
samples$samples[5]
samples$samples[8]

samples$samples[30]
samples$samples[31]

samples$samples[64]
samples$samples[67]

samples$samples[65]
samples$samples[66]

## Third order
samples$samples[11]
samples$samples[14]

samples$samples[38]
samples$samples[40]

ggplot(data, aes(x = a, y = b, color = KING, size = KING)) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Pairwise Relatedness", x = "Individual 1", y = "Individual 2", color = "Relatedness") +
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 1,size=15),
        axis.text.y=element_text(size=15),
        axis.title = element_text(size=16))+
   geom_label_repel(data = pairs_to_label, aes(label = names), size = 4, fill = "white", alpha = 0.8) 


pairs_to_label$type <-"first_order"
pairs_to_label$type[pairs_to_label$KING<0.17] <- "second_order"

p<- ggplot(data, aes(KING)) +
  geom_histogram(alpha = 0.7,bins=40) +
  labs(title = "", x = "Relatedness (King coefficient)", y = "", color = "Relatedness") +
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 1,size=15),
        axis.text.y=element_text(size=15),
        axis.title = element_text(size=16),
        legend.position = "none")+
  coord_flip()+
  geom_label_repel(
    data = pairs_to_label,
    aes(x = KING, y = 0, label = names,color=type),  # 'SampleID' or whichever label
    inherit.aes = FALSE,
    nudge_y = 120,
    angle = 90,       # rotate vertical
    vjust = 1,     # shift text a bit below the axis or above
    size = 4,
    direction="x",
    force=5,
    box.padding = 0.5,
    point.padding=0.3,
    segment.size=1,
    segment.linetype="dashed")+
  scale_color_manual(values=c("red","purple"))


p
ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/Supplemental_fig_relatedness_King.png", 
       plot = p, dpi = 300, width = 8, height = 10, units = "in")


