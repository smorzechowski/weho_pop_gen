# Classify WEHE contigs from HiFi based on coverage and heterozygosity
# Sophie MacRae Orzechowski
# August 2024

library(dplyr)
library(ggplot2)
library(stringr)

setwd("~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE assembly/findZX/Nesoptilotis_primary")


# get list of autosomal contigs to estimate median depth for each sex
# https://stackoverflow.com/questions/67192258/cant-read-txt-in-r-with-hash-tag-in-the-txt-file

lengths <- read.table("coverage/nleucotis_yamma_366917_f.p_ctg.filter.10000.fasta.fai",header=F,sep="\t",comment.char="")

matches <- read.table("synteny_lastal/T_guttatus/bestMatch.list",header=F,sep="\t",comment.char = "")

# Look at the first twenty autosomes from Zebra Finch
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_003957565.2/

autosome_list <- c("NC_044211.2","NC_044212.2","NC_044213.2","NC_044214.2",
                   "NC_044215.2","NC_044216.2","NC_044218.2","NC_044219.2",
                   "NC_044220.2","NC_044221.2","NC_044222.2","NC_044223.2",
                   "NC_044224.2","NC_044225.2","NC_044226.2","NC_044227.2",
                   "NC_044228.2","NC_044229.2","NC_044230.2","NC_044231.2")


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


# Curate a list of contigs that are homologous to Chr 5 in Zebra Finch

matches_chr5 <- matches[matches$V8=="NC_044217.2",]
matches_chr5$window_hits <- matches_chr5$V3 - matches_chr5$V2
matches_chr5_order <- matches_chr5[order(matches_chr5$V1),]
matches_chr5_order_nondup <- matches_chr5_order[!duplicated(matches_chr5_order$V1),]

matches_chr5_summary <- matches_chr5 %>% select(V1,window_hits) %>%
  group_by(V1) %>%
  summarize(windowsum = sum(window_hits,na.rm=TRUE))%>%
  as.data.frame()

matches_chr5_summary <- left_join(matches_chr5_summary,lengths)
matches_chr5_summary$percent_match <- matches_chr5_summary$windowsum/matches_chr5_summary$V2


# Designate a threshold for number of windows with hits or percent match
matches_chr5_final <- matches_chr5_summary[matches_chr5_summary$windowsum>35000,]
matches_chr5_final <- matches_chr5_summary[matches_chr5_summary$percent_match>0.10,]

# Curate a list of contigs that are homologous to Z in Zebra Finch

matches_Z <- matches[matches$V8=="NC_044241.2",]
matches_Z$window_hits <- matches_Z$V3 - matches_Z$V2
matches_Z_order <- matches_Z[order(matches_Z$V1),]
matches_Z_order_nondup <- matches_Z_order[!duplicated(matches_Z_order$V1),]

matches_Z_summary <- matches_Z %>% select(V1,window_hits) %>%
  group_by(V1) %>%
  summarize(windowsum = sum(window_hits,na.rm=TRUE))%>%
  as.data.frame()

matches_Z_summary <- left_join(matches_Z_summary,lengths)
matches_Z_summary$percent_match <- round(matches_Z_summary$windowsum/matches_Z_summary$V2,3)

# Designate a threshold for number of windows with hits or percent match
matches_Z_final <- matches_Z_summary[matches_Z_summary$windowsum>30000,]
#matches_Z_final <- matches_Z_summary[matches_Z_summary$percent_match>0.06,]

# Curate a list of contigs that are homologous to W in Zebra Finch

matches_W <- matches[matches$V8=="NC_045028.1",]
matches_W$window_hits <- matches_W$V3 - matches_W$V2
matches_W_order <- matches_W[order(matches_W$V1),]
matches_W_order_nondup <- matches_W_order[!duplicated(matches_W_order$V1),]


matches_W_summary <- matches_W %>% select(V1,window_hits) %>%
  group_by(V1) %>%
  summarize(windowsum = sum(window_hits,na.rm=TRUE))%>%
  as.data.frame()

matches_W_summary <- left_join(matches_W_summary,lengths)
matches_W_summary$percent_match <- round(matches_W_summary$windowsum/matches_W_summary$V2,3)

# Designate a threshold for number of windows with hits or percent match
matches_W_final <- matches_W_summary[matches_W_summary$windowsum>10000,]
#matches_W_final <- matches_W_summary[matches_W_summary$percent_match>0.06,]


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

#### Chr 5 coverage calculations ####

# Subset coverage based on window hits to Chr 5
coverage_chr5 <- coverage[coverage$V1%in%matches_chr5_final$V1,]
coverage_chr5$contig <- str_remove_all(coverage_chr5$V1,"nleucotis_yamma_366917_f#pri#ptg")


ggplot(coverage_chr5,aes(contig,log2FM))+
  geom_boxplot()


chr5_autosomal <- c("nleucotis_yamma_366917_f#pri#ptg000058l",
                    "nleucotis_yamma_366917_f#pri#ptg000069l",
                    "nleucotis_yamma_366917_f#pri#ptg000098l")

# In checking the matches, I've found that nleucotis_yamma_366917_f#pri#ptg000069l is also matching Chr 18 - NC_044230.2
# This might very well mean that a microchromosome fusion has occurred with truncated Chr 5. I will leave for now.

chr5_autosomal_lengths <- lengths[lengths$V1%in%chr5_autosomal,]
sum(chr5_autosomal_lengths$V2)

chr5_Z <- c("nleucotis_yamma_366917_f#pri#ptg000009l",
            "nleucotis_yamma_366917_f#pri#ptg000026l")


chr5_Z_lengths <- lengths[lengths$V1%in%chr5_Z,]
sum(chr5_Z_lengths$V2)

chr5_W <- c("nleucotis_yamma_366917_f#pri#ptg000039l",
            "nleucotis_yamma_366917_f#pri#ptg000071l")

chr5_W_lengths <- lengths[lengths$V1%in%chr5_W,]
sum(chr5_W_lengths$V2) 

# The result of 24605896 for W-linked Chr 5 means that much of this haplotype is missing from the primary assembly

#### Chr Z coverage calculations ####

# Subset coverage based on window hits to Chr Z
coverage_Z <- coverage[coverage$V1%in%matches_Z_final$V1,]
coverage_Z$contig <- str_remove_all(coverage_Z$V1,"nleucotis_yamma_366917_f#pri#ptg")

ggplot(coverage_Z,aes(contig,log2FM))+
  geom_boxplot()

chrZ_hits <- c("nleucotis_yamma_366917_f#pri#ptg000025l",
               "nleucotis_yamma_366917_f#pri#ptg000026l",
               "nleucotis_yamma_366917_f#pri#ptg000065l",
               "nleucotis_yamma_366917_f#pri#ptg000066l",
               "nleucotis_yamma_366917_f#pri#ptg000072l",
               "nleucotis_yamma_366917_f#pri#ptg000097l",
               "nleucotis_yamma_366917_f#pri#ptg000111l",
               "nleucotis_yamma_366917_f#pri#ptg000264l")

chrZ_hits_lengths <- lengths[lengths$V1%in%chrZ_hits,]
sum(chrZ_hits_lengths$V2)



chrW_hits <- c("nleucotis_yamma_366917_f#pri#ptg000030l",
               "nleucotis_yamma_366917_f#pri#ptg000039l",
               "nleucotis_yamma_366917_f#pri#ptg000062l",
               "nleucotis_yamma_366917_f#pri#ptg000071l",
               "nleucotis_yamma_366917_f#pri#ptg000088l",
               "nleucotis_yamma_366917_f#pri#ptg000213l")

chrW_hits_lengths <- lengths[lengths$V1%in%chrW_hits,]
sum(chrW_hits_lengths$V2)

#### Chr W coverage calculations ####

# Subset coverage based on window hits to Chr W
coverage_W <- coverage[coverage$V1%in%matches_W_final$V1,]
coverage_W$contig <- str_remove_all(coverage_W$V1,"nleucotis_yamma_366917_f#pri#ptg")

ggplot(coverage_W,aes(contig,log2FM))+
  geom_boxplot()+
  theme(axis.text.x=element_text(angle=45))


# adding to ChrW_hits from above

chrW_hits <- c("nleucotis_yamma_366917_f#pri#ptg000030l",
               "nleucotis_yamma_366917_f#pri#ptg000031l",
               "nleucotis_yamma_366917_f#pri#ptg000039l",
               "nleucotis_yamma_366917_f#pri#ptg000062l",
               "nleucotis_yamma_366917_f#pri#ptg000071l",
               "nleucotis_yamma_366917_f#pri#ptg000085l",
               "nleucotis_yamma_366917_f#pri#ptg000088l",
               "nleucotis_yamma_366917_f#pri#ptg000090l",
               "nleucotis_yamma_366917_f#pri#ptg000091l",
               "nleucotis_yamma_366917_f#pri#ptg000096l",
               "nleucotis_yamma_366917_f#pri#ptg000213l",
               "nleucotis_yamma_366917_f#pri#ptg001422l")

# ptg000030 is mostly W
# ptg000031 is mostly W
# ptg000039 is both W and 5 - already identified in Chr 5 above
# ptg000062 is mostly W
# ptg000071 is both W and 5 - already identified in Chr 5 above
# ptg000085 is mostly W and unplaced
# ptg000088 is mostly W
# ptg000090 is mostly W and unplaced
# ptg000091 is mostly unplaced and some W
# ptg000096 is mostly unplaced and some W
# ptg000213 is mostly W
# ptg001422 is W


# final coverage estimates for W
coverage_W_final <- coverage[coverage$V1%in%chrW_hits,]
coverage_W_final$contig <- str_remove_all(coverage_W_final$V1,"nleucotis_yamma_366917_f#pri#ptg")

ggplot(coverage_W_final,aes(contig,log2FM))+
  geom_boxplot()+
  theme(axis.text.x=element_text(angle=45))


chrW_hits_lengths <- lengths[lengths$V1%in%chrW_hits,]
sum(chrW_hits_lengths$V2) 

# The total amount of W-linked contigs I've found so far is 53636085, which means I'm still missing ~ 30 Mb of neo-W haplotype
# in this primary assembly

# adding to ChrZ hits from above
chrZ_hits <- c("nleucotis_yamma_366917_f#pri#ptg000025l",
               "nleucotis_yamma_366917_f#pri#ptg000026l",
               "nleucotis_yamma_366917_f#pri#ptg000057l",
               "nleucotis_yamma_366917_f#pri#ptg000065l",
               "nleucotis_yamma_366917_f#pri#ptg000066l",
               "nleucotis_yamma_366917_f#pri#ptg000072l",
               "nleucotis_yamma_366917_f#pri#ptg000097l",
               "nleucotis_yamma_366917_f#pri#ptg000111l",
               "nleucotis_yamma_366917_f#pri#ptg000107l",
               "nleucotis_yamma_366917_f#pri#ptg000264l")


chrZ_hits_lengths <- lengths[lengths$V1%in%chrZ_hits,]
sum(chrZ_hits_lengths$V2)

# final coverage estimates for Z
coverage_Z_final <- coverage[coverage$V1%in%chrZ_hits,]
coverage_Z_final$contig <- str_remove_all(coverage_Z_final$V1,"nleucotis_yamma_366917_f#pri#ptg")

ggplot(coverage_Z_final,aes(contig,log2FM))+
  geom_boxplot()+
  theme(axis.text.x=element_text(angle=45))


chrZ_hits_lengths <- lengths[lengths$V1%in%chrZ_hits,]
sum(chrZ_hits_lengths$V2) 


# look at lengths for Chr 5 Z-linked hits 
chr5_Z <- c("nleucotis_yamma_366917_f#pri#ptg000009l",
            "nleucotis_yamma_366917_f#pri#ptg000026l")


chr5_Z_lengths <- lengths[lengths$V1%in%chr5_Z,]
sum(chr5_Z_lengths$V2)

neoZ_hits <- c("nleucotis_yamma_366917_f#pri#ptg000025l",
                "nleucotis_yamma_366917_f#pri#ptg000026l",
                "nleucotis_yamma_366917_f#pri#ptg000057l",
                "nleucotis_yamma_366917_f#pri#ptg000065l",
                "nleucotis_yamma_366917_f#pri#ptg000066l",
                "nleucotis_yamma_366917_f#pri#ptg000072l",
                "nleucotis_yamma_366917_f#pri#ptg000097l",
                "nleucotis_yamma_366917_f#pri#ptg000111l",
                "nleucotis_yamma_366917_f#pri#ptg000107l",
                "nleucotis_yamma_366917_f#pri#ptg000264l",
                "nleucotis_yamma_366917_f#pri#ptg000009l")

neoZ_lengths <- lengths[lengths$V1%in%neoZ_hits,]
sum(neoZ_lengths$V2)

# ptg000025 is mostly ancestral Z
# ptg000026 is ancestral Z and Chr 5
# ptg000057 is mostly Chr 1 - why is log2FM so low? remove from final hits
# ptg000065 is ancestral Z
# ptg000066 is mostly ancestral Z
# ptg000072 is mostly ancestral Z
# ptg000097 is mostly ancestral Z
# ptg000111 is mostly ancestral Z
# ptg000107 is mostly W and Chr 35
# ptg000264 is mostly ancestral Z
# ptg000009 is Chr 5

neoZ_hits_final <- c("nleucotis_yamma_366917_f#pri#ptg000025l",
                     "nleucotis_yamma_366917_f#pri#ptg000026l",
                     "nleucotis_yamma_366917_f#pri#ptg000517l", # 517 is a scaffold identified by ragtag
                     "nleucotis_yamma_366917_f#pri#ptg000065l",
                     "nleucotis_yamma_366917_f#pri#ptg000066l",
                     "nleucotis_yamma_366917_f#pri#ptg000072l",
                     "nleucotis_yamma_366917_f#pri#ptg000097l",
                     "nleucotis_yamma_366917_f#pri#ptg000111l",
                     "nleucotis_yamma_366917_f#pri#ptg000264l",
                     "nleucotis_yamma_366917_f#pri#ptg000009l")

neoZ_lengths <- lengths[lengths$V1%in%neoZ_hits_final,]
sum(neoZ_lengths$V2)

#### Calculate mean depth and log2FM depth #####

coverage_summary <- coverage %>% select(V1,Fnorm,Mnorm,V4,V5) %>%
  group_by(V1) %>%
  summarize(meanFnorm = mean(Fnorm,na.rm=TRUE),
            meanMnorm = mean(Mnorm,na.rm=TRUE),
            meanF = mean(V4,na.rm=TRUE),
            meanM = mean(V5,na.rm=TRUE))%>%
  as.data.frame()


coverage_summary$log2FM <- log2(coverage_summary$meanF/coverage_summary$meanM)

### Additional primary contigs found that are likely W-linked and have no synteny matches ####
# nleucotis_yamma_366917_f#pri#ptg001094l
# nleucotis_yamma_366917_f#pri#ptg001195l
# nleucotis_yamma_366917_f#pri#ptg000318l
# nleucotis_yamma_366917_f#pri#ptg001220l
# 

