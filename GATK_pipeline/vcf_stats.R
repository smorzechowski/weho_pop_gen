# VCFtools summary of unfiltered vcf files

# S M Orzechowski
# 9 November 2024


library(ggplot2)
library(forcats)
library(stringr)
library(colorBlindness)
library(RColorBrewer)
library(dplyr)
library(tidyr)

setwd("~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen")

##### Individual depth #####
depth <- read.table("Nleu_interval_1_combined_varsite.idepth",header=T)
depth$sex <- str_remove(depth$INDV,'[:alpha:]+_[0-9]+_')
depth$sex[depth$sex=='U']='M'

ggplot(depth,aes(fct_reorder(INDV,MEAN_DEPTH),MEAN_DEPTH,fill=sex))+
  geom_bar(stat='identity')+
  xlab('Sample')+
  ylab('Mean depth')+
  theme(axis.text.x=element_blank(),
        axis.title=element_text(size=16),
        axis.text.y=element_text(size=14))+
  scale_fill_brewer(palette="Paired")+
  geom_text(aes(label=INDV), vjust = 0,angle=90)
  #scale_fill_manual()

ggplot(depth,aes(sex,MEAN_DEPTH,fill=sex))+
  geom_violin()+
  ylab('Mean depth')+
  xlab('Sex')+
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=14))+
  scale_fill_brewer(palette="Paired")


##### Phred quality scores ######
qual <- read.table("Nleu_interval_1_combined_varsite.lqual",header=T)

qual_subset <- qual[sample(nrow(qual),40000),]

ggplot(qual_subset[qual_subset$QUAL<500,],aes(QUAL))+
  geom_density(fill = "dodgerblue1",color="black",alpha=0.3)+
  geom_vline(xintercept=30)+
  xlab('Quality (Phred score)')+
  ylab('Density')+
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=14))

##### Variant-level depth ######
var_depth <- read.table("Nleu_interval_1_combined_varsite.ldepth.mean",header=T)
var_subset <- var_depth[sample(nrow(var_depth),50000),]


ggplot(var_subset[var_subset$VAR_DEPTH<100,],aes(VAR_DEPTH))+
  geom_density(fill = "dodgerblue1",color="black",alpha=0.3)+
  geom_vline(xintercept=20)+
  geom_vline(xintercept=11,linetype=2)+
  geom_vline(xintercept=3)+
  xlab('Variant depth')+
  ylab('Density')+
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=14))


##### Missingness per variant site ######

var_miss <- read.table("Nleu_interval_1_combined_varsite.lmiss",header=T)
var_miss_subset <- var_miss[sample(nrow(var_miss),500000),]

summary(var_miss_subset$F_MISS)

ggplot(var_miss_subset,aes(F_MISS))+
  geom_density(fill = "dodgerblue1",color="black",alpha=0.3)+
  xlab('Genotype missingness (% individuals)')+
  geom_vline(xintercept=0.1,linetype=2)+
  ylab('Density')+
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=14))


##### Missingness per individual #####

ind_miss  <- read.table("Nleu_interval_1_combined_varsite.imiss", sep = "\t",
                        col.names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)

ggplot(ind_miss,aes(fct_reorder(ind,fmiss),fmiss))+
  geom_bar(stat='identity')+
  xlab('Sample')+
  ylab('Proportion missing data')+
  theme(axis.text.x=element_blank(),
        axis.title=element_text(size=16),
        axis.text.y=element_text(size=14))+
  scale_fill_brewer(palette="Paired")
#scale_fill_manual()


##### Variant frequency ######
var_freq <- read.table("Nleu_interval_1_combined_varsite.frq",col.names=c('chr',"pos","nalleles","nchr","a1","a2"),skip=1)
var_freq_subset <- var_freq[sample(nrow(var_freq),500000),]

# find minor allele frequency
var_freq_subset$maf <- var_freq_subset %>% select(a1, a2) %>% apply(1, function(z) min(z))

ggplot(var_freq_subset,aes(maf))+
  geom_density(fill = "dodgerblue1",color="black",alpha=0.3)+
  geom_vline(xintercept=0.05,linetype=2)+
  xlab('Minor allele frequency')+
  ylab('Density')+
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=14))

###### Relatedness
relate <- read.table("interval_1.relatedness2",header=T)
relate <- relate[,c(1,2,7)]

relate_wide <- pivot_wider(relate,names_from=INDV1,values_from=RELATEDNESS_PHI)
relate_wide <- data.frame(relate_wide)
row.names(relate_wide) <- relate_wide$INDV2
relate_wide <- relate_wide[,-1]

data <- data.matrix(relate_wide)
heatmap(data)


