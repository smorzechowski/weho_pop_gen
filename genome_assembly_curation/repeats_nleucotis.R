# RepeatMasker summary N. leucotis
# S Orzechowski

library(ggplot2)
library(dplyr)


##### interval information
#scaffold_5_RagTag 0 26782956 # added W
#scaffold_5_RagTag 26782957 32320175 # ancestral W
#scaffold_5_RagTag 32320176 55423182 # added W
#scaffold_5_RagTag 55423183 73135152 # ancestral W  
#scaffold_5_RagTag 73135153 90514351 # new PAR -- WHICH SHOULD BE HARDMASKED!!!!!= 
#newPAR<-90514351 - 73135153
#scaffold_2 total 129027458 
#scaffold_2 New PAR boundary 129027458 - 17379198 = 111648260
# Scaffold_2 ancestral Z boundary: 74301111 ? or thereabouts or now...
#####

setwd("~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE assembly/Nleu_final_genome_softmask")

data <- read.csv('Nleucotis_hifi_v1.0_hc_sm_fx_scaffolded_PAR_masked.fasta.out.csv',header=F)


neoZ <- data[data$V5=='scaffold_2_RagTag',]
neoW <- data[data$V5=='scaffold_5_RagTag',]

neoW$region[neoW$V7<26782956] ='addedW'
neoW$region[neoW$V7>26782956 & neoW$V7<32320175] ='ancestralW'
neoW$region[neoW$V7>32320175 & neoW$V7<55423182] ='addedW'
neoW$region[neoW$V7>55423183 & neoW$V7<73135152] ='ancestralW'
neoW$reptot <- neoW$V7-neoW$V6

Wtotal <- neoW %>% group_by(V12,region)%>%
  summarise(total_repeat=sum(reptot))%>%
  data.frame()

ancestralW_1 <- 32320175 - 26782956
ancestralW_2 <- 73135152 - 55423183
ancestralW_total <- ancestralW_1 + ancestralW_2 

addedW_1 <- 26782956
addedW_2 <- 55423182 - 32320175

addedW_total <- addedW_1 + addedW_2

Wtotal$total_length=NA
Wtotal$total_length[Wtotal$region=="ancestralW"] <- 23249188
Wtotal$total_length[Wtotal$region=="addedW"] = 49885963



## neoZ
neoZ$region[neoZ$V7>111648260] ='newPAR'
neoZ$region[neoZ$V7>74301111 & neoZ$V7<111648260] ='addedZ'
neoZ$region[neoZ$V7<74301111] ='ancestralZ'
neoZ$reptot <- neoZ$V7-neoZ$V6
Ztotal <- neoZ %>% group_by(V12,region)%>%
  summarise(total_repeat=sum(reptot))%>%
  data.frame()

Ztotal$total_length <- 0
Ztotal$total_length[Ztotal$region=="newPAR"] <- 17379198
addedZ = 111648260 - 74301111
Ztotal$total_length[Ztotal$region=="addedZ"] = 37347149
Ztotal$total_length[Ztotal$region=="ancestralZ"] = 74301111


#Ztotal$percent <- round(Ztotal$total_repeat/Ztotal$total_length,4)

#data$Elements<- factor(data$Elements,levels=c("DNA","LINEs","LTRs","Low_complexity","Simple_repeats","Unclassified"))

combined <- rbind(Wtotal,Ztotal)
combined$region <- factor(combined$region,levels=c('ancestralW','addedW','ancestralZ','addedZ','newPAR'))
colnames(combined)[1] <- "Class"

#ggplot(data,aes(reorder(factor(Region),Order),Masked_length/Total_length,fill=Elements))+
ggplot(combined,aes(factor(region),total_repeat/total_length,fill=Class))+
  geom_bar(stat="identity",position='stack')+
  # scale_color_viridis(discrete=TRUE)+
  #scale_fill_viridis(discrete=TRUE)+
  xlab("")+
  #xlab("Region")+
  ylim(0,1)+
  #ylab("Proportion of TEs per region")+
  ylab("")+
  theme_classic()+
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=15),
        axis.text.x =element_text(size=15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=20))


