# Classify WEHE contigs from HiFi based on coverage and heterozygosity
# Sophie MacRae Orzechowski

# August 2024

library(dplyr)
library(ggplot2)
library(stringr)

setwd("~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE assembly/findZX/Nesoptilotis_hap2")


# get list of autosomal contigs to estimate median depth for each sex
# https://stackoverflow.com/questions/67192258/cant-read-txt-in-r-with-hash-tag-in-the-txt-file

lengths <- read.table("coverage/nleucotis_yamma_366917_f.hap2.p_ctg.filter.10000.fasta.fai",header=F,sep="\t",comment.char="")

matches <- read.table("synteny_lastal/E_cyanotis/bestMatch.list",header=F,sep="\t",comment.char = "")

# Look at the first twenty autosomes from Blue-faced Honeyeater


autosome_list <- c("scaffold_1","scaffold_3","scaffold_4","scaffold_6",
                   "scaffold_7","scaffold_8","scaffold_9","scaffold_10",
                   "scaffold_11","scaffold_12","scaffold_13","scaffold_14",
                   "scaffold_15","scaffold_16","scaffold_17","scaffold_18",
                   "scaffold_19","scaffold_20","scaffold_21","scaffold_22")


# Curate a list of contigs that are homologous to these autosomes in Zebra Finch

matches_auto <- matches[matches$V8%in%autosome_list,]
matches_auto$window_hits <- matches_auto$V3 - matches_auto$V2
matches_auto_order <- matches_auto[order(matches_auto$V1),]
matches_auto_order_nondup <- matches_auto_order[!duplicated(matches_auto_order$V1),]


matches_auto_summary <- matches_auto %>% select(V1,window_hits) %>%
  group_by(V1) %>%
  summarize(windowsum = sum(window_hits,na.rm=TRUE))%>%
  as.data.frame()

matches_auto_summary <- left_join(matches_auto_summary,lengths)


# Curate a list of contigs that are homologous to neo-Z in Blue-faced Honeyeater

matches_neoZ <- matches[matches$V8=="scaffold_2",]
matches_neoZ$window_hits <- matches_neoZ$V3 - matches_neoZ$V2
matches_neoZ_order <- matches_neoZ[order(matches_neoZ$V1),]
matches_neoZ_order_nondup <- matches_neoZ_order[!duplicated(matches_neoZ_order$V1),]

matches_neoZ_summary <- matches_neoZ %>% select(V1,window_hits) %>%
  group_by(V1) %>%
  summarize(windowsum = sum(window_hits,na.rm=TRUE))%>%
  as.data.frame()

matches_neoZ_summary <- left_join(matches_neoZ_summary,lengths)
matches_neoZ_summary$percent_match <- matches_neoZ_summary$windowsum/matches_neoZ_summary$V2
# Designate a threshold for number of windows with hits, or percentage match
matches_neoZ_final <- matches_neoZ_summary[matches_neoZ_summary$windowsum>35000,]
matches_neoZ_final <- matches_neoZ_summary[matches_neoZ_summary$percent_match>0.10,]


# Curate a list of contigs that are homologous to neo-W in Blue-faced Honeyeater

matches_neoW <- matches[matches$V8=="scaffold_5",]
matches_neoW$window_hits <- matches_neoW$V3 - matches_neoW$V2
matches_neoW_order <- matches_neoW[order(matches_neoW$V1),]
matches_neoW_order_nondup <- matches_neoW_order[!duplicated(matches_neoW_order$V1),]


matches_neoW_summary <- matches_neoW %>% select(V1,window_hits) %>%
  group_by(V1) %>%
  summarize(windowsum = sum(window_hits,na.rm=TRUE))%>%
  as.data.frame()

matches_neoW_summary <- left_join(matches_neoW_summary,lengths)
matches_neoW_summary$percent_match <- matches_neoW_summary$windowsum/matches_neoW_summary$V2
# Designate a threshold for number of windows with hits or percent match
matches_neoW_final <- matches_neoW_summary[matches_neoW_summary$windowsum>10000,]
matches_neoW_final <- matches_neoW_summary[matches_neoW_summary$percent_match>0.10,]


##### Estimate coverage for each genomic compartment


coverage <- read.table("coverage/gencov.mismatch.0.0.out",header=F,sep="\t",comment.char = "")

# Calculate median depth for selected 20 autosomes for each sex

auto_coverage <- coverage[coverage$V1%in%matches_auto_order_nondup$V1,]

# females
median_sampleV4 <- median(auto_coverage$V4)

# males
median_sampleV5 <- median(auto_coverage$V5)

# calculate log2FM following Schield et al. 2022

coverage$Fnorm <- coverage$V4/median_sampleV4
coverage$Mnorm <- coverage$V5/median_sampleV5

coverage$log2FM <- log2(coverage$Fnorm/coverage$Mnorm)

#### neo-Z coverage calculations ####

# Subset coverage based on window hits to Chr 5
coverage_neoZ <- coverage[coverage$V1%in%matches_neoZ_final$V1,]
coverage_neoZ$contig <- str_remove_all(coverage_neoZ$V1,"nleucotis_yamma_366917_f.hap2#pri#")


ggplot(coverage_neoZ,aes(contig,log2FM))+
  geom_boxplot()


neoZ <- c("nleucotis_yamma_366917_f.hap2#pri#h2tg000008l",
          "nleucotis_yamma_366917_f.hap2#pri#h2tg000050l",
          "nleucotis_yamma_366917_f.hap2#pri#h2tg000076l",
          "nleucotis_yamma_366917_f.hap2#pri#h2tg000082l",
          "nleucotis_yamma_366917_f.hap2#pri#h2tg000223l",
          "nleucotis_yamma_366917_f.hap2#pri#h2tg000289l",
          "nleucotis_yamma_366917_f.hap2#pri#h2tg000772l")

# "nleucotis_yamma_366917_f.hap2#pri#h2tg000023l", this is a repeat, I've discovered!! Very similar contig in hap1 and hap2


neoZ_lengths <- lengths[lengths$V1%in%neoZ,]
neoZ_hap2 <- sum(neoZ_lengths$V2)
neoZ_hap2


#### neo-W coverage calculations ####

# Subset coverage based on window hits to Chr W
coverage_neoW <- coverage[coverage$V1%in%matches_neoW_final$V1,]

coverage_neoW$contig <- str_remove_all(coverage_neoW$V1,"nleucotis_yamma_366917_f.hap2#pri#")

ggplot(coverage_neoW,aes(contig,log2FM))+
  geom_boxplot()+
  theme(axis.text.x=element_text(angle=45))

chrW_hits <- c("nleucotis_yamma_366917_f.hap2#pri#h2tg000039l",
               "nleucotis_yamma_366917_f.hap2#pri#h2tg000069l",
               "nleucotis_yamma_366917_f.hap2#pri#h2tg000081l",
               "nleucotis_yamma_366917_f.hap2#pri#h2tg000101l",
               "nleucotis_yamma_366917_f.hap2#pri#h2tg000287l",
               "nleucotis_yamma_366917_f.hap2#pri#h2tg000769l",
               "nleucotis_yamma_366917_f.hap2#pri#h2tg001519l",
               "nleucotis_yamma_366917_f.hap2#pri#h2tg001715l")

chrW_hits_lengths <- lengths[lengths$V1%in%chrW_hits,]
neoW_hap2 <- sum(chrW_hits_lengths$V2)
neoW_hap2


#### Calculate mean depth and log2FM depth #####

summary <- coverage %>% select(V1,Fnorm,Mnorm,V4,V5) %>%
  group_by(V1) %>%
  summarize(meanFnorm = mean(Fnorm,na.rm=TRUE),
            meanMnorm = mean(Mnorm,na.rm=TRUE),
            meanF = mean(V4,na.rm=TRUE),
            meanM = mean(V5,na.rm=TRUE))%>%
  as.data.frame()


summary$log2FM <- log2(summary$meanF/summary$meanM)
