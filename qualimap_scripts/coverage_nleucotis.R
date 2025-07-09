# Assign sex to Nesoptilotis leucotis samples based on coverage
# Sophie MacRae Orzechowski
# December 2024


library(ggplot2)

setwd("~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE assembly")


# all_coverage_results_cleaned_final.txt contains the two individuals from Moonbi sequenced seperately

data <- read.table("all_coverage_results_cleaned_final.txt",header=F,sep=" ",comment.char = "")

data_filt <- data[data$V1=="scaffold_5_RagTag" |
                  data$V1=="scaffold_2_RagTag" |
                  data$V1=="scaffold_1_RagTag",]


ggplot(data_filt,aes(V6,V4))+
  geom_bar(stat='identity')+
  facet_grid(.~V1)


data_W <- data[data$V1=="scaffold_5_RagTag",]
ggplot(data_W,aes(V6,V4))+
  geom_bar(stat='identity')+
  coord_flip()+
  ylab('Neo-W Coverage (excluding new PAR)')+
  xlab("")


data_1a <- data[data$V1=="scaffold_1_RagTag",]
ggplot(data_1a,aes(V6,V4))+
  geom_bar(stat='identity')+
  coord_flip()+
  ylab('1A coverage (autosomal)')+
  xlab("")
