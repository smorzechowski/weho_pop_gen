# GEA hits

#
library(ggplot2)

setwd("~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/LEA")


data <- read.csv('GEA_hits_by_chromosome.csv',header=T)

data$Chr <- factor(data$Chromosome,levels=c("Chr_10","Ancestral","Chr_13","ptg15l","Chr_4",
                                            "Added","Chr_12","Chr_17","ptg45l","Chr_11",
                                            "Chr_19","Chr_3","Chr_8","Chr_29","Chr_9","Chr_1","Chr_14",
                                            "Chr_15","Chr_18","Chr_20","Chr_21","Chr_26","Chr_6","Chr_7"))

data$Category <- factor(data$Category,levels=c("Ancestral Z","Added Z","Autosomes"))

p <- ggplot(data,aes(Chr,GEA,fill=Category))+
  geom_bar(stat='identity',color='black')+
  ylab('# GEA hits')+
  xlab('Chromosome')+
  theme_bw()+
  theme(axis.text.x=element_text(angle=60,size=13),
        axis.text.y=element_text(size=13),
        axis.title=element_text(size=14),
        legend.title=element_blank(),
        legend.position='top',
        legend.text=element_text(size=14))+
  scale_fill_manual(values=c("black","white","grey"))


p 
ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/GEA_hits_chromosome.png", 
       plot = p, dpi = 350, width = 8, height = 6, units = "in")