# Make synteny plots of structural rearrangements in WEHE
# Sophie MacRae Orzechowski
# February 2025
###################################################################################

library(RIdeogram)
library(pafr)
library(stringr)
library(Cairo)
library(svglite)
library(rsvg)
library(dplyr)

setwd("~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE assembly/paf")
#source("ideogram_abs_length_updated2.R")
#source("ideogram.R")
#data(human_karyotype, package="RIdeogram")
#head(human_karyotype)
#data(karyotype_dual_comparison, package="RIdeogram")
#head(karyotype_dual_comparison)
#data(synteny_dual_comparison, package="RIdeogram")
#head(synteny_dual_comparison)
###################################################################################


paf_df <- read_paf("scaffold_10_22_RagTag_scaffold_10_22.paf")

# Create a data frame for the query information
output_q <- paf_df %>%
  mutate(
    Chr = qname,    # Use query name
    Start = 1,      # Start is fixed at 1
    End = qlen,     # End is the query length
    fill = "969696",      # Empty string for fill
    species = "WEHE",
    size = 12,      # Set size to 12
    color = "252525"      # Empty string for color
  ) %>%
  select(Chr, Start, End, fill, species, size, color)

# Create a data frame for the target information
output_t <- paf_df %>%
  mutate(
    Chr = tname,    # Use target name
    Start = 1,      # Start is fixed at 1
    End = tlen,     # End is the target length
    fill = "636363",      # Empty string for fill
    species ="BFHE",
    size = 12,      # Set size to 12
    color = "252525"      # Empty string for color
  ) %>%
  select(Chr, Start, End, fill, species, size, color)

# Combine the two data frames into one
karyotype_df <- bind_rows(output_t, output_q) %>%
  distinct()%>%
  as.data.frame()



paf_df$qname <- as.numeric(str_remove_all(paf_df$qname,"scaffold_|_RagTag"))
paf_df$tname <- as.numeric(str_remove_all(paf_df$tname,"scaffold_|_RagTag"))

paf_df$qname <- as.numeric(as.factor(paf_df$qname))
paf_df$tname <- as.numeric(as.factor(paf_df$tname))


# Create a data frame for synteny
synteny_df <- paf_df %>%
  mutate(
    Species_1 = tname,    # Use target name
    Start_1 = tstart,      # Start is fixed at 1
    End_1 = tend,     # End is the target length
    Species_2 = qname,      # Empty string for fill
    Start_2 =qstart,
    End_2 =qend,      # Set size to 12
    fill = "cccccc"      # Empty string for color
  ) %>%
  select(Species_1, Start_1, End_1,Species_2,Start_2,End_2,fill)%>%
  as.data.frame()

# color where the rearrangement is 
synteny_df$fill[synteny_df$Start_1>9675047] <- "0000FF"


ideogram(karyotype = karyotype_df, synteny = synteny_df,output="scaffold_10_22_RagTag.svg")
convertSVG("scaffold_10_22_RagTag.svg", file="scaffold_10_22_RagTag", device = "png")



###################################################################################


paf_df <- read_paf("prig000047l_scaffold_8_12_15_21.paf")

# Create a data frame for the query information
output_q <- paf_df %>%
  mutate(
    Chr = qname,    # Use query name
    Start = 1,      # Start is fixed at 1
    End = qlen,     # End is the query length
    fill = "636363",      # Empty string for fill 00FF00
    species = "WEHE",
    size = 12,      # Set size to 12
    color = "252525"      # Empty string for color
  ) %>%
  select(Chr, Start, End, fill, species, size, color)

# Create a data frame for the target information
output_t <- paf_df %>%
  mutate(
    Chr = tname,    # Use target name
    Start = 1,      # Start is fixed at 1
    End = tlen,     # End is the target length
    fill = "636363",      # Empty string for fill
    species ="BFHE",
    size = 12,      # Set size to 12
    color = "252525"      # Empty string for color
  ) %>%
  select(Chr, Start, End, fill, species, size, color)

# Combine the two data frames into one
karyotype_df <- bind_rows(output_t, output_q) %>%
  distinct()%>%
  as.data.frame()

#karyotype_df$fill[karyotype_df$Chr=="scaffold_15"] <- "00FF00"
#karyotype_df$fill[karyotype_df$Chr=="scaffold_8"] <- "FFA500"


#paf_df$qname <- as.numeric(str_remove_all(paf_df$qname,"scaffold_|_RagTag"))
#paf_df$tname <- as.numeric(str_remove_all(paf_df$tname,"scaffold_|_RagTag"))

paf_df$qname <- as.numeric(as.factor(paf_df$qname))


karyotype_df <- karyotype_df[c(3,2,1,4,5),]
paf_df$tname[paf_df$tname=="scaffold_8"] <- 1
paf_df$tname[paf_df$tname=="scaffold_12"] <- 2
paf_df$tname[paf_df$tname=="scaffold_15"] <- 3
paf_df$tname[paf_df$tname=="scaffold_21"] <- 4 
paf_df$tname <- as.numeric(paf_df$tname)

# Create a data frame for synteny
synteny_df <- paf_df %>%
  mutate(
    Species_1 = tname,    # Use target name
    Start_1 = tstart,      # Start is fixed at 1
    End_1 = tend,     # End is the target length
    Species_2 = qname,      # Empty string for fill
    Start_2 =qstart,
    End_2 =qend,      # Set size to 12
    fill = "cccccc"      # Empty string for color
  ) %>%
  select(Species_1, Start_1, End_1,Species_2,Start_2,End_2,fill)%>%
  as.data.frame()



#synteny_df$Start_1[synteny_df$Species_1==1] <- 38036437 - synteny_df$Start_1[synteny_df$Species_1==1] 
#synteny_df$End_1[synteny_df$Species_1==1] <- 38036437 - synteny_df$End_1[synteny_df$Species_1==1] 


#synteny_df$Start_new[synteny_df$Species_1==1] <- synteny_df$End_1[synteny_df$Species_1==1]
#synteny_df$End_new[synteny_df$Species_1==1] <- synteny_df$Start_1[synteny_df$Species_1==1]

#synteny_df$Start_1[synteny_df$Species_1==1] <- synteny_df$Start_new[synteny_df$Species_1==1]
#synteny_df$End_1[synteny_df$Species_1==1] <- synteny_df$End_new[synteny_df$Species_1==1]


#synteny_df <- select(synteny_df, -Start_new, -End_new)

#start_test = 37874992 -15100000
#end_test = 55974910 - 15100000

#synteny_df$fill[synteny_df$Species_1==1] <- "0000FF"
#synteny_df$fill[synteny_df$Start_2>] <- "0000FG"
#synteny_df$fill[synteny_df$Start_2>start_test & synteny_df$Start_2<end_test] <- "0000FG"


ideogram(karyotype = karyotype_df, synteny = synteny_df,output="contig_047.svg")
convertSVG("contig_047.svg", file="contig_047", device = "png")


##########################################################################################
# Look at hap1 contig that corroborates the fusion between 12 and 15 in primary contig 047!!!! (involving Chr_8)

paf_df <- read_paf("h1tg_000098l_scaff_15_12_names_updated.paf")

# Create a data frame for the query information
output_q <- paf_df %>%
  mutate(
    Chr = qname,    # Use query name
    Start = 1,      # Start is fixed at 1
    End = qlen,     # End is the query length
    fill = "636363",      # Empty string for fill 00FF00
    species = "WEHE",
    size = 12,      # Set size to 12
    color = "252525"      # Empty string for color
  ) %>%
  select(Chr, Start, End, fill, species, size, color)

# Create a data frame for the target information
output_t <- paf_df %>%
  mutate(
    Chr = tname,    # Use target name
    Start = 1,      # Start is fixed at 1
    End = tlen,     # End is the target length
    fill = "636363",      # Empty string for fill
    species ="BFHE",
    size = 12,      # Set size to 12
    color = "252525"      # Empty string for color
  ) %>%
  select(Chr, Start, End, fill, species, size, color)

# Combine the two data frames into one
karyotype_df <- bind_rows(output_t, output_q) %>%
  distinct()%>%
  as.data.frame()


paf_df$qname <- as.numeric(as.factor(paf_df$qname))


paf_df$tname[paf_df$tname=="scaffold_12"] <- 1
paf_df$tname[paf_df$tname=="scaffold_15"] <- 2

paf_df$tname <- as.numeric(paf_df$tname)

# Create a data frame for synteny
synteny_df <- paf_df %>%
  mutate(
    Species_1 = tname,    # Use target name
    Start_1 = tstart,      # Start is fixed at 1
    End_1 = tend,     # End is the target length
    Species_2 = qname,      # Empty string for fill
    Start_2 =qstart,
    End_2 =qend,      # Set size to 12
    fill = "cccccc"      # Empty string for color
  ) %>%
  select(Species_1, Start_1, End_1,Species_2,Start_2,End_2,fill)%>%
  as.data.frame()


ideogram(karyotype = karyotype_df, synteny = synteny_df,output="h1tg_000098l_scaff_15_12.svg")
convertSVG("h1tg_000098l_scaff_15_12.svg", file="h1tg_000098l_scaff_15_12", device = "png")



#############################################################################################
# look at fusion between scaffold_8 and 21 and 12!

paf_df <- read_paf("h1tg000054_scaff_8_21_12_names_updated.paf")

# Create a data frame for the query information
output_q <- paf_df %>%
  mutate(
    Chr = qname,    # Use query name
    Start = 1,      # Start is fixed at 1
    End = qlen,     # End is the query length
    fill = "636363",      # Empty string for fill 00FF00
    species = "WEHE",
    size = 12,      # Set size to 12
    color = "252525"      # Empty string for color
  ) %>%
  select(Chr, Start, End, fill, species, size, color)

# Create a data frame for the target information
output_t <- paf_df %>%
  mutate(
    Chr = tname,    # Use target name
    Start = 1,      # Start is fixed at 1
    End = tlen,     # End is the target length
    fill = "636363",      # Empty string for fill
    species ="BFHE",
    size = 12,      # Set size to 12
    color = "252525"      # Empty string for color
  ) %>%
  select(Chr, Start, End, fill, species, size, color)

# Combine the two data frames into one
karyotype_df <- bind_rows(output_t, output_q) %>%
  distinct()%>%
  as.data.frame()


paf_df$qname <- as.numeric(as.factor(paf_df$qname))

paf_df$tname[paf_df$tname=="scaffold_21"] <- 1
paf_df$tname[paf_df$tname=="scaffold_8"] <- 2
paf_df$tname[paf_df$tname=="scaffold_12"] <- 3

paf_df$tname <- as.numeric(paf_df$tname)

# Create a data frame for synteny
synteny_df <- paf_df %>%
  mutate(
    Species_1 = tname,    # Use target name
    Start_1 = tstart,      # Start is fixed at 1
    End_1 = tend,     # End is the target length
    Species_2 = qname,      # Empty string for fill
    Start_2 =qstart,
    End_2 =qend,      # Set size to 12
    fill = "cccccc"      # Empty string for color
  ) %>%
  select(Species_1, Start_1, End_1,Species_2,Start_2,End_2,fill)%>%
  as.data.frame()


ideogram(karyotype = karyotype_df, synteny = synteny_df,output="h1tg000054_scaff_8_21_12.svg")
convertSVG("h1tg000054_scaff_8_21_12.svg", file="h1tg000054_scaff_8_21_12", device = "png")







###################################################################################
# look at potential inversion on HAP2 contig 000014l

paf_df <- read_paf("hap2_000014l_scaffold_9_names_updated.paf")


# Compute midpoints for query and target alignment blocks
paf_df$query_mid  <- (paf_df$qstart + paf_df$qend) / 2
paf_df$target_mid <- (paf_df$tstart + paf_df$tend) / 2

# Some measure of similarity - need to check on this
paf_df$percentID = paf_df$nmatch / paf_df$alen
queryStartTemp = paf_df$qstart

# Flip starts, ends for negative strand alignments
paf_df$qstart[which(paf_df$strand == "-")] = paf_df$qend[which(paf_df$strand == "-")]
paf_df$qend[which(paf_df$strand == "-")] = queryStartTemp[which(paf_df$strand == "-")]
rm(queryStartTemp)

#cat(paste0("\nNumber of alignments: ", nrow(alignments),"\n"))
#cat(paste0("Number of query sequences: ", length(unique(alignments$queryID)),"\n"))



paf_df_filt <- paf_df[-c(13:16),]

# Create a basic dot plot with ggplot2
ggplot(paf_df,aes(x=qstart/1e6,y=tstart/1e6)) +
  geom_segment(data=paf_df,aes(x=tstart/1e6,xend=tend/1e6,y=qstart/1e6,yend=qend/1e6),linewidth = 2)+
  #  geom_point(alpha = 0.7) +
  ylim(0,17.7)+
  xlim(0,36)+
  #geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  labs(y = "Query position (Mbp)", x = "Target position (Mbp)") +
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15))






# Create a basic dot plot with ggplot2
#ggplot(paf_df_filt, aes(x = query_mid/1e6, y = target_mid/1e6, color = strand)) +
#  geom_segment(data=paf_df_filt,aes(x=qstart/1e6,xend=qend/1e6,y=tstart/1e6,yend=tend/1e6))+
#  geom_point(alpha = 0.7) +
#  xlim(0,17.7)+
#  geom_line()+
#  #geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
#  labs(x = "Query position (bp)", y = "Target position (bp)", color = "Strand") +
#  theme_minimal()


# Create a data frame for the query information
output_q <- paf_df %>%
  mutate(
    Chr = qname,    # Use query name
    Start = 1,      # Start is fixed at 1
    End = qlen,     # End is the query length
    fill = "969696",      # Empty string for fill
    species = "WEHE",
    size = 12,      # Set size to 12
    color = "252525"      # Empty string for color
  ) %>%
  select(Chr, Start, End, fill, species, size, color)

# Create a data frame for the target information
output_t <- paf_df %>%
  mutate(
    Chr = tname,    # Use target name
    Start = 1,      # Start is fixed at 1
    End = tlen,     # End is the target length
    fill = "636363",      # Empty string for fill
    species ="BFHE",
    size = 12,      # Set size to 12
    color = "252525"      # Empty string for color
  ) %>%
  select(Chr, Start, End, fill, species, size, color)

# Combine the two data frames into one
karyotype_df <- bind_rows(output_t, output_q) %>%
  distinct()%>%
  as.data.frame()


paf_df$qname <- as.numeric(as.factor(paf_df$qname))
paf_df$tname <- as.numeric(as.factor(paf_df$tname))

# Create a data frame for synteny
synteny_df <- paf_df %>%
  mutate(
    Species_1 = tname,    # Use target name
    Start_1 = tstart,      # Start is fixed at 1
    End_1 = tend,     # End is the target length
    Species_2 = qname,      # Empty string for fill
    Start_2 =qstart,
    End_2 =qend,      # Set size to 12
    fill = "cccccc"      # Empty string for color
  ) %>%
  select(Species_1, Start_1, End_1,Species_2,Start_2,End_2,fill)%>%
  as.data.frame()

synteny_df$fill[synteny_df$Start_2<10864590] <- "0000FF"


ideogram(karyotype = karyotype_df, synteny = synteny_df,output="hap2_014_inversion.svg")
convertSVG("hap2_014_inversion.svg", file="hap2_014_inversion", device = "png")





###################################################################################
# look at potential inversion on contig priptg000015l #REFERENCE# REVERSE COMPLEMENTED

paf_df <- read_paf("query15_against_Chr9_RC.paf")

# Create a data frame for the query information
output_q <- paf_df %>%
  mutate(
    Chr = qname,    # Use query name
    Start = 1,      # Start is fixed at 1
    End = qlen,     # End is the query length
    fill = "969696",      # Empty string for fill
    species = "WEHE",
    size = 12,      # Set size to 12
    color = "252525"      # Empty string for color
  ) %>%
  select(Chr, Start, End, fill, species, size, color)

# Create a data frame for the target information
output_t <- paf_df %>%
  mutate(
    Chr = tname,    # Use target name
    Start = 1,      # Start is fixed at 1
    End = tlen,     # End is the target length
    fill = "636363",      # Empty string for fill
    species ="BFHE",
    size = 12,      # Set size to 12
    color = "252525"      # Empty string for color
  ) %>%
  select(Chr, Start, End, fill, species, size, color)

# Combine the two data frames into one
karyotype_df <- bind_rows(output_t, output_q) %>%
  distinct()%>%
  as.data.frame()


paf_df$qname <- as.numeric(as.factor(paf_df$qname))
paf_df$tname <- as.numeric(as.factor(paf_df$tname))


# Create a data frame for synteny
synteny_df <- paf_df %>%
  mutate(
    Species_1 = tname,    # Use target name
    Start_1 = tstart,      # Start is fixed at 1
    End_1 = tend,     # End is the target length
    Species_2 = qname,      # Empty string for fill
    Start_2 =qstart,
    End_2 =qend,      # Set size to 12
    fill = "cccccc"      # Empty string for color
  ) %>%
  select(Species_1, Start_1, End_1,Species_2,Start_2,End_2,fill)%>%
  as.data.frame()

synteny_df$fill[synteny_df$Start_2<11044568] <- "0000FF"


ideogram(karyotype = karyotype_df, synteny = synteny_df,output="contig_015_REF_RC_inversion.svg")
convertSVG("contig_015_RC_inversion.svg", file="contig_015_REF_RC_inversion", device = "png")




###################################################################################
# look at potential inversion on contig priptg000015l REVERSE COMPLEMENTED

paf_df <- read_paf("query15_against_ref9.paf")

# Some measure of similarity - need to check on this
paf_df$percentID = paf_df$nmatch / paf_df$alen
queryStartTemp = paf_df$qstart

# Flip starts, ends for negative strand alignments
paf_df$qstart[which(paf_df$strand == "-")] = paf_df$qend[which(paf_df$strand == "-")]
paf_df$qend[which(paf_df$strand == "-")] = queryStartTemp[which(paf_df$strand == "-")]
rm(queryStartTemp)



# Create a data frame for the query information
output_q <- paf_df %>%
  mutate(
    Chr = qname,    # Use query name
    Start = 1,      # Start is fixed at 1
    End = qlen,     # End is the query length
    fill = "969696",      # Empty string for fill
    species = "WEHE",
    size = 12,      # Set size to 12
    color = "252525"      # Empty string for color
  ) %>%
  select(Chr, Start, End, fill, species, size, color)

# Create a data frame for the target information
output_t <- paf_df %>%
  mutate(
    Chr = tname,    # Use target name
    Start = 1,      # Start is fixed at 1
    End = tlen,     # End is the target length
    fill = "636363",      # Empty string for fill
    species ="BFHE",
    size = 12,      # Set size to 12
    color = "252525"      # Empty string for color
  ) %>%
  select(Chr, Start, End, fill, species, size, color)

# Combine the two data frames into one
karyotype_df <- bind_rows(output_t, output_q) %>%
  distinct()%>%
  as.data.frame()


paf_df$qname <- as.numeric(as.factor(paf_df$qname))
paf_df$tname <- as.numeric(as.factor(paf_df$tname))


# Create a data frame for synteny
synteny_df <- paf_df %>%
  mutate(
    Species_1 = tname,    # Use target name
    Start_1 = tstart,      # Start is fixed at 1
    End_1 = tend,     # End is the target length
    Species_2 = qname,      # Empty string for fill
    Start_2 =qstart,
    End_2 =qend,      # Set size to 12
    fill = "cccccc"      # Empty string for color
  ) %>%
  select(Species_1, Start_1, End_1,Species_2,Start_2,End_2,fill)%>%
  as.data.frame()

synteny_df$fill[synteny_df$Start_2<11044568] <- "0000FF"


ideogram(karyotype = karyotype_df, synteny = synteny_df,output="contig_015_RC_inversion2.svg")
convertSVG("contig_015_RC_inversion2.svg", file="contig_015_RC_inversion2", device = "png")



###################################################################################
# look at potential inversion on contig priptg000015l!!!!!!

paf_df <- read_paf("scaffold_9_015l.paf")


# Compute midpoints for query and target alignment blocks
paf_df$query_mid  <- (paf_df$qstart + paf_df$qend) / 2
paf_df$target_mid <- (paf_df$tstart + paf_df$tend) / 2

paf_df_filt <- paf_df[-c(17,19,28:32,34,36),]




# Some measure of similarity - need to check on this
paf_df$percentID = paf_df$nmatch / paf_df$alen
queryStartTemp = paf_df$qstart

# Flip starts, ends for negative strand alignments
paf_df$qstart[which(paf_df$strand == "-")] = paf_df$qend[which(paf_df$strand == "-")]
paf_df$qend[which(paf_df$strand == "-")] = queryStartTemp[which(paf_df$strand == "-")]
rm(queryStartTemp)


#paf_df_filt <- paf_df[-c(13:16),]

# Create a basic dot plot with ggplot2
p <- ggplot(paf_df,aes(x=qstart/1e6,y=tstart/1e6)) +
  geom_segment(data=paf_df,aes(x=tstart/1e6,xend=tend/1e6,y=qstart/1e6,yend=qend/1e6),linewidth = 2)+
  #  geom_point(alpha = 0.7) +
  ylim(0,17.7)+
  xlim(0,36)+
  #geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  labs(y = "Query position (Mbp)", x = "Target position (Mbp)") +
  theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=15))+
  scale_x_reverse()

p

ggsave("C:/Users/sophi/Documents/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/figures/priptg000015l_dotplot_inversion.png", 
       plot = p, dpi = 350, width = 6, height = 6, units = "in")



####################################################################################
# Create a basic dot plot with ggplot2
#ggplot(paf_df_filt, aes(x = query_mid/1e6, y = target_mid/1e6, color = strand)) +
#  geom_point(alpha = 0.7) +
#  geom_segment(data=paf_df_filt,aes(x=qstart/1e6,xend=qend/1e6,y=tstart/1e6,yend=tend/1e6))+
#  ylim(0,36)+
#  xlim(0,17.7)+
#  geom_line()+
# # geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
#  labs(x = "Query position (bp)", y = "Target position (bp)", color = "Strand") +
#  theme_minimal()


#ggplot(paf_df_filt, aes(x = query_mid/1e6, y = target_mid/1e6, color = strand)) +
#  geom_point(alpha = 0.7) +
#  geom_segment(data=paf_df_filt,aes(x=qstart/1e6,xend=qend/1e6,y=tstart/1e6,yend=tend/1e6))+
#  ylim(0,36)+
#  xlim(0,17.7)+
# # geom_line()+
#  # geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
#  labs(x = "Query position (bp)", y = "Target position (bp)", color = "Strand") +
#  theme_minimal()



# Create a data frame for the query information
output_q <- paf_df %>%
  mutate(
    Chr = qname,    # Use query name
    Start = 1,      # Start is fixed at 1
    End = qlen,     # End is the query length
    fill = "969696",      # Empty string for fill
    species = "WEHE",
    size = 12,      # Set size to 12
    color = "252525"      # Empty string for color
  ) %>%
  select(Chr, Start, End, fill, species, size, color)

# Create a data frame for the target information
output_t <- paf_df %>%
  mutate(
    Chr = tname,    # Use target name
    Start = 1,      # Start is fixed at 1
    End = tlen,     # End is the target length
    fill = "636363",      # Empty string for fill
    species ="BFHE",
    size = 12,      # Set size to 12
    color = "252525"      # Empty string for color
  ) %>%
  select(Chr, Start, End, fill, species, size, color)

# Combine the two data frames into one
karyotype_df <- bind_rows(output_t, output_q) %>%
  distinct()%>%
  as.data.frame()

#karyotype_df$Start[karyotype_df$Chr=="priptg000015l"] <- 17909547
#karyotype_df$End[karyotype_df$Chr=="priptg000015l"] <- 1

#paf_df$qname <- as.numeric(str_remove_all(paf_df$qname,"scaffold_|_RagTag"))
#paf_df$tname <- as.numeric(str_remove_all(paf_df$tname,"scaffold_|_RagTag"))

paf_df$qname <- as.numeric(as.factor(paf_df$qname))
paf_df$tname <- as.numeric(as.factor(paf_df$tname))


# 35259631 -
# 17909547 - 
# Create a data frame for synteny
synteny_df <- paf_df %>%
  mutate(
    Species_1 = tname,    # Use target name
    Start_1 = 35259631 - tstart,      # Start is fixed at 1
    End_1 = 35259631 - tend,     # End is the target length
    Species_2 = qname,      # Empty string for fill
    Start_2 =qstart,
    End_2 =qend,      # Set size to 12
    fill = "cccccc"      # Empty string for color
  ) %>%
  select(Species_1, Start_1, End_1,Species_2,Start_2,End_2,fill)%>%
  as.data.frame()

synteny_df$fill[synteny_df$Start_2>6664979] <- "0000FF"
#synteny_df$Start_1[synteny_df$Start_2>6964979] <- "0000FF"

#synteny_df$Start_new <- synteny_df$End_2
#synteny_df$End_new <- synteny_df$Start_2

#synteny_df$Start_2 <- synteny_df$Start_new
#synteny_df$End_2 <- synteny_df$End_new
#

#synteny_df <- select(synteny_df, -Start_new, -End_new)



ideogram(karyotype = karyotype_df, synteny = synteny_df,output="contig_015_inversion_flippedtest.svg")
convertSVG("contig_015_inversion_flippedtest.svg", file="contig_015_inversion_flippedtest", device = "png")


# Create a data frame for synteny
synteny_df <- paf_df %>%
  mutate(
    Species_1 = tname,    # Use target name
    Start_1 = tstart,      # Start is fixed at 1
    End_1 = tend,     # End is the target length
    Species_2 = qname,      # Empty string for fill
    Start_2 =qstart,
    End_2 =qend,      # Set size to 12
    fill = "cccccc"      # Empty string for color
  ) %>%
  select(Species_1, Start_1, End_1,Species_2,Start_2,End_2,fill)%>%
  as.data.frame()

synteny_df$fill[synteny_df$Start_2>6664979] <- "0000FF"
#synteny_df$Start_1[synteny_df$Start_2>6964979] <- "0000FF"


ideogram(karyotype = karyotype_df, synteny = synteny_df,output="contig_015_inversion.svg")
convertSVG("contig_015_inversion.svg", file="contig_015_inversion", device = "png")


##############################################################################################
# Create a data frame for synteny
# NEED to manually flip the chromosome here!!!!!!!!
##############################################################################################

# NOTE:
# I wasn't sure if this is appropriate to reverse complement just a portion of the contig
# This is based on the results of manually reverse complementing and then creating the alignment again
# But it might actually be after all! See my conversation with o1 in Evernote


synteny_df <- paf_df %>%
  mutate(
    Species_1 = tname,    # Use target name
    Start_1 =  35259631-tstart,      # Start is fixed at 1
    End_1 = 35259631- tend,     # End is the target length
    Species_2 = qname,      # Empty string for fill
    Start_2 =qstart,
    End_2 = qend,      # Set size to 12
    fill = "cccccc"      # Empty string for color
  ) %>%
  select(Species_1, Start_1, End_1,Species_2,Start_2,End_2,fill)%>%
  as.data.frame()

synteny_df$fill[synteny_df$Start_2>6664979] <- "0000FF"


ideogram(karyotype = karyotype_df, synteny = synteny_df,output="contig_015_inversion_flipped_target.svg")
convertSVG("contig_015_inversion_flipped_target.svg", file="contig_015_inversion_flipped_target", device = "png")




### Flip target and reverse alignments 

synteny_df <- paf_df %>%
  mutate(
    Species_1 = tname,    # Use target name
    Start_1 =  35259631-tstart,      # Start is fixed at 1
    End_1 = 35259631- tend,     # End is the target length
    Species_2 = qname,      # Empty string for fill
    Start_2 =qstart,
    End_2 = qend,      # Set size to 12
    fill = "cccccc"      # Empty string for color
  ) %>%
  select(Species_1, Start_1, End_1,Species_2,Start_2,End_2,fill)%>%
  as.data.frame()

synteny_df$fill[synteny_df$Start_2>6664979] <- "0000FF"
synteny_df$Start_new[synteny_df$Start_2<6664979] <- synteny_df$End_1[synteny_df$Start_2<6664979]
synteny_df$End_new[synteny_df$Start_2<6664979] <- synteny_df$Start_1[synteny_df$Start_2<6664979]

synteny_df$Start_1[synteny_df$Start_2<6664979] <- synteny_df$Start_new[synteny_df$Start_2<6664979]
synteny_df$End_1[synteny_df$Start_2<6664979] <- synteny_df$End_new[synteny_df$Start_2<6664979]

synteny_df <- select(synteny_df, -Start_new, -End_new)


ideogram(karyotype = karyotype_df, synteny = synteny_df,output="contig_015_inversion_flipped_target_rev_align.svg")
convertSVG("contig_015_inversion_flipped_target_rev_align.svg", file="contig_015_inversion_flipped_target_rev_align", device = "png")





###################################################################################
# look at potential fusion!!!!!!

paf_df <- read_paf("scaffold_10_11_pri45l.paf")

# Create a data frame for the query information
output_q <- paf_df %>%
  mutate(
    Chr = qname,    # Use query name
    Start = 1,      # Start is fixed at 1
    End = qlen,     # End is the query length
    fill = "969696",      # Empty string for fill
    species = "WEHE",
    size = 12,      # Set size to 12
    color = "252525"      # Empty string for color
  ) %>%
  select(Chr, Start, End, fill, species, size, color)

# Create a data frame for the target information
output_t <- paf_df %>%
  mutate(
    Chr = tname,    # Use target name
    Start = 1,      # Start is fixed at 1
    End = tlen,     # End is the target length
    fill = "636363",      # Empty string for fill
    species ="BFHE",
    size = 12,      # Set size to 12
    color = "252525"      # Empty string for color
  ) %>%
  select(Chr, Start, End, fill, species, size, color)

# Combine the two data frames into one
karyotype_df <- bind_rows(output_t, output_q) %>%
  distinct()%>%
  as.data.frame()

karyotype_df <- karyotype_df[c(2,1,3),]

#karyotype_df$Start[karyotype_df$Chr=="priptg000015l"] <- 17909547
#karyotype_df$End[karyotype_df$Chr=="priptg000015l"] <- 1

#paf_df$qname <- as.numeric(str_remove_all(paf_df$qname,"scaffold_|_RagTag"))
#paf_df$tname <- as.numeric(str_remove_all(paf_df$tname,"scaffold_|_RagTag"))

paf_df$qname <- as.numeric(as.factor(paf_df$qname))
paf_df$tname <- as.numeric(as.factor(paf_df$tname))

paf_df$tname[paf_df$tname=="scaffold_10"] <- 1
paf_df$tname[paf_df$tname=="scaffold_11"] <- 2


# OK, this is not appropriate here, or at least it needs to be paired with reorienting RC as forward 
# if the majority of alignments are on the negative strand. This would erroneously display inversions if you did not do the 
# second step, which is what dotplotly does -- see the documentation. However, in the potential inversion, it doesn't really matter I think.
# The inversion signature will remain either way. This goes way deeper than I probably needed to, ah well.


# Some measure of similarity - need to check on this
#paf_df$percentID = paf_df$nmatch / paf_df$alen
#queryStartTemp = paf_df$qstart

# Flip starts, ends for negative strand alignments
#paf_df$qstart[which(paf_df$strand == "-")] = paf_df$qend[which(paf_df$strand == "-")]
#paf_df$qend[which(paf_df$strand == "-")] = queryStartTemp[which(paf_df$strand == "-")]
#rm(queryStartTemp)


# Create a data frame for synteny
synteny_df <- paf_df %>%
  mutate(
    Species_1 = tname,    # Use target name
    Start_1 = tstart,      # Start is fixed at 1
    End_1 = tend,     # End is the target length
    Species_2 =  qname,      # Empty string for fill
    Start_2 =  qstart,
    End_2 =qend,      # Set size to 12
    fill = "cccccc"      # Empty string for color
  ) %>%
  select(Species_1, Start_1, End_1,Species_2,Start_2,End_2,fill)%>%
  as.data.frame()

synteny_df$fill[synteny_df$Start_2<12834860] <- "0000FF"
synteny_df$Start_1[synteny_df$Species_1==1] <- 30612747 - synteny_df$Start_1[synteny_df$Species_1==1] 
synteny_df$End_1[synteny_df$Species_1==1] <- 30612747 - synteny_df$End_1[synteny_df$Species_1==1] 

synteny_df$Start_1[synteny_df$Species_1==2] <- 25293088 - synteny_df$Start_1[synteny_df$Species_1==2] 
synteny_df$End_1[synteny_df$Species_1==2] <- 25293088 - synteny_df$End_1[synteny_df$Species_1==2] 


synteny_df$Start_new <- synteny_df$End_1
synteny_df$End_new <- synteny_df$Start_1

synteny_df$Start_1 <- synteny_df$Start_new
synteny_df$End_1 <- synteny_df$End_new


synteny_df <- select(synteny_df, -Start_new, -End_new)



ideogram(karyotype = karyotype_df, synteny = synteny_df,output="contig_045_fusion_11_reversed.svg")
convertSVG("contig_045_fusion_11_reversed.svg", file="contig_045_fusion_11_reversed", device = "png")


###################################################################################
# look at potential fusion!!!!!!

paf_df <- read_paf("scaffold_8_RagTag_8_12_15_21.paf")

# Create a data frame for the query information
output_q <- paf_df %>%
  mutate(
    Chr = qname,    # Use query name
    Start = 1,      # Start is fixed at 1
    End = qlen,     # End is the query length
    fill = "969696",      # Empty string for fill
    species = "WEHE",
    size = 12,      # Set size to 12
    color = "252525"      # Empty string for color
  ) %>%
  select(Chr, Start, End, fill, species, size, color)

# Create a data frame for the target information
output_t <- paf_df %>%
  mutate(
    Chr = tname,    # Use target name
    Start = 1,      # Start is fixed at 1
    End = tlen,     # End is the target length
    fill = "636363",      # Empty string for fill
    species ="BFHE",
    size = 12,      # Set size to 12
    color = "252525"      # Empty string for color
  ) %>%
  select(Chr, Start, End, fill, species, size, color)

# Combine the two data frames into one
karyotype_df <- bind_rows(output_t, output_q) %>%
  distinct()%>%
  as.data.frame()

paf_df$qname <- as.numeric(as.factor(paf_df$qname))

paf_df$tname[paf_df$tname=="scaffold_8"] <- 1
paf_df$tname[paf_df$tname=="scaffold_21"] <- 2 
paf_df$tname[paf_df$tname=="scaffold_12"] <- 3
paf_df$tname[paf_df$tname=="scaffold_15"] <- 4

paf_df$tname <- as.numeric(paf_df$tname)


# Create a data frame for synteny
synteny_df <- paf_df %>%
  mutate(
    Species_1 = tname,    # Use target name
    Start_1 = tstart,      # Start is fixed at 1
    End_1 = tend,     # End is the target length
    Species_2 = qname,      # Empty string for fill
    Start_2 =qstart,
    End_2 =qend,      # Set size to 12
    fill = "cccccc"      # Empty string for color
  ) %>%
  select(Species_1, Start_1, End_1,Species_2,Start_2,End_2,fill)%>%
  as.data.frame()

#synteny_df$Start_1[synteny_df$Species_1==1] <- 	38036437 - synteny_df$Start_1[synteny_df$Species_1==1] 
#synteny_df$End_1[synteny_df$Species_1==1] <- 	38036437 - synteny_df$End_1[synteny_df$Species_1==1] 

#scaffold_8_RagTag 1 10724985
#scaffold_8_RagTag 37874992 56196169

synteny_df$fill[synteny_df$Start_2<10724985] <- "0000FF"
synteny_df$fill[synteny_df$Start_2>37874992& synteny_df$Start_2<56196169] <- "0000FF"

#synteny_df$Start_new[synteny_df$Species_1==1] <- synteny_df$End_1[synteny_df$Species_1==1]
#synteny_df$End_new[synteny_df$Species_1==1] <- synteny_df$Start_1[synteny_df$Species_1==1]

#synteny_df$Start_1[synteny_df$Species_1==1] <- synteny_df$Start_new[synteny_df$Species_1==1]
#synteny_df$End_1[synteny_df$Species_1==1] <- synteny_df$End_new[synteny_df$Species_1==1]



#synteny_df$Start_new <- synteny_df$End_2
#synteny_df$End_new <- synteny_df$Start_2

#synteny_df$Start_2 <- synteny_df$Start_new
#synteny_df$End_2 <- synteny_df$End_new


#synteny_df <- select(synteny_df, -Start_new, -End_new)


# Flip scaffold 12 so that we can see that there's no inversion -- just oriented on the negative strand
synteny_df$Start_new[synteny_df$Species_1==3] <- 25079814 - synteny_df$End_1[synteny_df$Species_1==3]
synteny_df$End_new[synteny_df$Species_1==3] <- 25079814 - synteny_df$Start_1[synteny_df$Species_1==3]

synteny_df$Start_1[synteny_df$Species_1==3] <- synteny_df$Start_new[synteny_df$Species_1==3]
synteny_df$End_1[synteny_df$Species_1==3] <- synteny_df$End_new[synteny_df$Species_1==3]



# Don't bother flipping scaffold_15 because it's in the same orientation in the other figures

synteny_df <- select(synteny_df, -Start_new, -End_new)



ideogram(karyotype = karyotype_df, synteny = synteny_df,output="scaffold_8_RagTag_fusion_test.svg")
convertSVG("scaffold_8_RagTag_fusion_test.svg", file="scaffold_8_RagTag_fusion_test", device = "png")


###################################################################################
# look at potential fusion!!!!!!

paf_df <- read_paf("prtg000003l_scaffold_13_15.paf")

# Create a data frame for the query information
output_q <- paf_df %>%
  mutate(
    Chr = qname,    # Use query name
    Start = 1,      # Start is fixed at 1
    End = qlen,     # End is the query length
    fill = "969696",      # Empty string for fill
    species = "WEHE",
    size = 12,      # Set size to 12
    color = "252525"      # Empty string for color
  ) %>%
  select(Chr, Start, End, fill, species, size, color)

# Create a data frame for the target information
output_t <- paf_df %>%
  mutate(
    Chr = tname,    # Use target name
    Start = 1,      # Start is fixed at 1
    End = tlen,     # End is the target length
    fill = "636363",      # Empty string for fill
    species ="BFHE",
    size = 12,      # Set size to 12
    color = "252525"      # Empty string for color
  ) %>%
  select(Chr, Start, End, fill, species, size, color)

# Combine the two data frames into one
karyotype_df <- bind_rows(output_t, output_q) %>%
  distinct()%>%
  as.data.frame()


karyotype_df <- karyotype_df[c(2,1,3),]

#karyotype_df$Start[karyotype_df$Chr=="priptg000015l"] <- 17909547
#karyotype_df$End[karyotype_df$Chr=="priptg000015l"] <- 1

#paf_df$qname <- as.numeric(str_remove_all(paf_df$qname,"scaffold_|_RagTag"))
#paf_df$tname <- as.numeric(str_remove_all(paf_df$tname,"scaffold_|_RagTag"))

paf_df$qname <- as.numeric(as.factor(paf_df$qname))
#paf_df$tname <- as.numeric(as.factor(paf_df$tname))

paf_df$tname[paf_df$tname=="scaffold_13"] <- 1
paf_df$tname[paf_df$tname=="scaffold_15"] <- 2

paf_df$tname <- as.numeric(paf_df$tname)


#20080491 - 
#20080491 - 
# Create a data frame for synteny
synteny_df <- paf_df %>%
  mutate(
    Species_1 = tname,    # Use target name
    Start_1 = tstart,      # Start is fixed at 1
    End_1 = tend,     # End is the target length
    Species_2 = qname,      # Empty string for fill
    Start_2 = 20080491 - qstart,
    End_2 = 20080491 - qend,      # Set size to 12
    fill = "cccccc"      # Empty string for color
  ) %>%
  select(Species_1, Start_1, End_1,Species_2,Start_2,End_2,fill)%>%
  as.data.frame()

synteny_df$fill[synteny_df$Start_2>12834860] <- "0000FF"
#synteny_df$Start_new[synteny_df$Species_1==2] <- 20886537 - synteny_df$End_1[synteny_df$Species_1==2]
#synteny_df$End_new[synteny_df$Species_1==2] <- 20886537 - synteny_df$Start_1[synteny_df$Species_1==2]

#synteny_df$Start_1[synteny_df$Species_1==2] <- synteny_df$Start_new[synteny_df$Species_1==2]
#synteny_df$End_1[synteny_df$Species_1==2] <- synteny_df$End_new[synteny_df$Species_1==2]

#synteny_df <- select(synteny_df, -Start_new, -End_new)


synteny_df$Start_new <- synteny_df$End_2
synteny_df$End_new <- synteny_df$Start_2

synteny_df$Start_2 <- synteny_df$Start_new
synteny_df$End_2 <- synteny_df$End_new


synteny_df <- select(synteny_df, -Start_new, -End_new)



ideogram(karyotype = karyotype_df, synteny = synteny_df,output="prtg00003l_scaffold_13_15_fusion.svg")
convertSVG("prtg00003l_scaffold_13_15_fusion.svg", file="prtg00003l_scaffold_13_15_fusion", device = "png")


###################################################################################
# look at potential fusion!!!!!!

paf_df <- read_paf("scaffold_15_RagTag_scaffold_15_11_37.paf")


# Create a data frame for the query information
output_q <- paf_df %>%
  mutate(
    Chr = qname,    # Use query name
    Start = 1,      # Start is fixed at 1
    End = qlen,     # End is the query length
    fill = "969696",      # Empty string for fill
    species = "WEHE",
    size = 12,      # Set size to 12
    color = "252525"      # Empty string for color
  ) %>%
  select(Chr, Start, End, fill, species, size, color)

# Create a data frame for the target information
output_t <- paf_df %>%
  mutate(
    Chr = tname,    # Use target name
    Start = 1,      # Start is fixed at 1
    End = tlen,     # End is the target length
    fill = "636363",      # Empty string for fill
    species ="BFHE",
    size = 12,      # Set size to 12
    color = "252525"      # Empty string for color
  ) %>%
  select(Chr, Start, End, fill, species, size, color)

# Combine the two data frames into one
karyotype_df <- bind_rows(output_t, output_q) %>%
  distinct()%>%
  as.data.frame()


karyotype_df <- karyotype_df[c(2,1,3,4),]

paf_df$qname <- as.numeric(as.factor(paf_df$qname))
paf_df$tname <- as.numeric(as.factor(paf_df$tname))

#paf_df$tname[paf_df$tname=="scaffold_15"] <- 1
#paf_df$tname[paf_df$tname=="scaffold_11"] <- 2


# Create a data frame for synteny
synteny_df <- paf_df %>%
  mutate(
    Species_1 = tname,    # Use target name
    Start_1 = tstart,      # Start is fixed at 1
    End_1 = tend,     # End is the target length
    Species_2 = qname,      # Empty string for fill
    Start_2 = qstart,
    End_2 = qend,      # Set size to 12
    fill = "cccccc"      # Empty string for color
  ) %>%
  select(Species_1, Start_1, End_1,Species_2,Start_2,End_2,fill)%>%
  as.data.frame()


synteny_df$fill[synteny_df$Start_2<8774970] <- "0000FF"



ideogram(karyotype = karyotype_df, synteny = synteny_df,output="scaffold_15_RagTag.svg")
convertSVG("scaffold_15_RagTag.svg", file="scaffold_15_RagTag", device = "png")



###################################################################################
# look at potential fusion!!!!!!

paf_df <- read_paf("scaffold_6_RagTag_6_12.paf")


# Create a data frame for the query information
output_q <- paf_df %>%
  mutate(
    Chr = qname,    # Use query name
    Start = 1,      # Start is fixed at 1
    End = qlen,     # End is the query length
    fill = "969696",      # Empty string for fill
    species = "WEHE",
    size = 12,      # Set size to 12
    color = "252525"      # Empty string for color
  ) %>%
  select(Chr, Start, End, fill, species, size, color)

# Create a data frame for the target information
output_t <- paf_df %>%
  mutate(
    Chr = tname,    # Use target name
    Start = 1,      # Start is fixed at 1
    End = tlen,     # End is the target length
    fill = "636363",      # Empty string for fill
    species ="BFHE",
    size = 12,      # Set size to 12
    color = "252525"      # Empty string for color
  ) %>%
  select(Chr, Start, End, fill, species, size, color)

# Combine the two data frames into one
karyotype_df <- bind_rows(output_t, output_q) %>%
  distinct()%>%
  as.data.frame()


#karyotype_df <- karyotype_df[c(2,1,3,4),]

paf_df$qname <- as.numeric(as.factor(paf_df$qname))
#paf_df$tname <- as.numeric(as.factor(paf_df$tname))

paf_df$tname[paf_df$tname=="scaffold_6"] <- 1
paf_df$tname[paf_df$tname=="scaffold_12"] <- 2
paf_df$tname <- as.numeric(paf_df$tname)

# Create a data frame for synteny
synteny_df <- paf_df %>%
  mutate(
    Species_1 = tname,    # Use target name
    Start_1 = tstart,      # Start is fixed at 1
    End_1 = tend,     # End is the target length
    Species_2 = qname,      # Empty string for fill
    Start_2 = qstart,
    End_2 = qend,      # Set size to 12
    fill = "cccccc"      # Empty string for color
  ) %>%
  select(Species_1, Start_1, End_1,Species_2,Start_2,End_2,fill)%>%
  as.data.frame()


synteny_df$fill[synteny_df$Start_2>8125001&synteny_df$Start_2<32874990] <- "0000FF"
#scaffold_6_RagTag 8125001 32874990


ideogram(karyotype = karyotype_df, synteny = synteny_df,output="scaffold_6_RagTag.svg")
convertSVG("scaffold_6_RagTag.svg", file="scaffold_6_RagTag", device = "png")



###################################################################################
# Extra code
###################################################################################
# pafr example

plot_synteny(data,
             q_chrom = "scaffold_10_RagTag",t_chrom ="scaffold_10",centre=TRUE)

plot_synteny(data,
             q_chrom = "scaffold_10_RagTag",t_chrom ="scaffold_22",centre=TRUE)


# ideogram example
data(karyotype_dual_comparison, package="RIdeogram")
head(karyotype_dual_comparison)

data(synteny_dual_comparison, package="RIdeogram")
head(synteny_dual_comparison)
str(synteny_dual_comparison)
ideogram(karyotype = karyotype_dual_comparison, synteny = synteny_dual_comparison)
convertSVG("chromosome.svg", device = "png")


# convert svg with different packages

# Open a Cairo-based SVG device:
CairoSVG(filename = "chromosome_test.svg", width = 8.2677, height = 11.6929)
# Create your plot here. For example:
plot(1:10, 1:10, main = "Example Plot")
# When done, close the device:
dev.off()

# Open the svglite device:
svglite("chromosome_test.svg", width = 8.2677, height = 11.6929)
# Create your plot:
plot(1:10, 1:10, main = "Example Plot")
# Close the device:
dev.off()
rsvg_png("chromosome_test.svg", "chromosome_testing.png", width = 8.2677 * 300, height = 11.6929 * 300)


# try using geom_curve (did not work)
ggplot(data) +
  geom_curve(aes(x = tstart, y = qstart, 
                 xend = tend, yend = qend),
             curvature = 0.2, 
             color = "blue", 
             size = 1) +
  labs(x = "Reference position", y = "Query position")





