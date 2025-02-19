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


ideogram(karyotype = karyotype_df, synteny = synteny_df,output="chromosome_test.svg")
convertSVG("chromosome_test.svg", device = "png")



###################################################################################


paf_df <- read_paf("prig000047l_scaffold_8_12_15_21.paf")

# Create a data frame for the query information
output_q <- paf_df %>%
  mutate(
    Chr = qname,    # Use query name
    Start = 1,      # Start is fixed at 1
    End = qlen,     # End is the query length
    fill = "00FF00",      # Empty string for fill
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

karyotype_df$fill[karyotype_df$Chr=="scaffold_15"] <- "00FF00"
karyotype_df$fill[karyotype_df$Chr=="scaffold_8"] <- "FFA500"


#paf_df$qname <- as.numeric(str_remove_all(paf_df$qname,"scaffold_|_RagTag"))
#paf_df$tname <- as.numeric(str_remove_all(paf_df$tname,"scaffold_|_RagTag"))

paf_df$qname <- as.numeric(as.factor(paf_df$qname))
#paf_df$tname <- as.numeric(as.factor(paf_df$tname))

paf_df$tname[paf_df$tname=="scaffold_15"] <- 1
paf_df$tname[paf_df$tname=="scaffold_12"] <- 2
paf_df$tname[paf_df$tname=="scaffold_8"] <- 3
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



synteny_df$Start_1[synteny_df$Species_1==3] <- 38036437 - synteny_df$Start_1[synteny_df$Species_1==3] 
synteny_df$End_1[synteny_df$Species_1==3] <- 38036437 - synteny_df$End_1[synteny_df$Species_1==3] 


synteny_df$Start_new[synteny_df$Species_1==3] <- synteny_df$End_1[synteny_df$Species_1==3]
synteny_df$End_new[synteny_df$Species_1==3] <- synteny_df$Start_1[synteny_df$Species_1==3]

synteny_df$Start_1[synteny_df$Species_1==3] <- synteny_df$Start_new[synteny_df$Species_1==3]
synteny_df$End_1[synteny_df$Species_1==3] <- synteny_df$End_new[synteny_df$Species_1==3]


synteny_df <- select(synteny_df, -Start_new, -End_new)



synteny_df$fill[synteny_df$Species_1==3] <- "0000FF"


ideogram(karyotype = karyotype_df, synteny = synteny_df,output="chromosome_047.svg")
convertSVG("chromosome_047.svg", file="chromosome_047", device = "png")



###################################################################################
# look at potential inversion!!!!!!

paf_df <- read_paf("scaffold_9_015l.paf")

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


# Create a data frame for synteny
synteny_df <- paf_df %>%
  mutate(
    Species_1 = tname,    # Use target name
    Start_1 = 35259631 - tend,      # Start is fixed at 1
    End_1 = 35259631 - tstart,     # End is the target length
    Species_2 = qname,      # Empty string for fill
    Start_2 =qstart,
    End_2 =qend,      # Set size to 12
    fill = "cccccc"      # Empty string for color
  ) %>%
  select(Species_1, Start_1, End_1,Species_2,Start_2,End_2,fill)%>%
  as.data.frame()

synteny_df$fill[synteny_df$Start_2>6664979] <- "0000FF"
#synteny_df$Start_1[synteny_df$Start_2>6964979] <- "0000FF"


ideogram(karyotype = karyotype_df, synteny = synteny_df,output="chromosome_015_inversion_flipped.svg")
convertSVG("chromosome_015_inversion_flipped.svg", file="chromosome_015_inversion_flipped", device = "png")


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


ideogram(karyotype = karyotype_df, synteny = synteny_df,output="chromosome_015_inversion.svg")
convertSVG("chromosome_015_inversion.svg", file="chromosome_015_inversion", device = "png")


##############################################################################################
# Create a data frame for synteny
# NEED to manually flip the chromosome here!!!!!!!!
##############################################################################################

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
synteny_df$Start_new[synteny_df$Start_2<6664979] <- synteny_df$End_1[synteny_df$Start_2<6664979]
synteny_df$End_new[synteny_df$Start_2<6664979] <- synteny_df$Start_1[synteny_df$Start_2<6664979]

synteny_df$Start_1[synteny_df$Start_2<6664979] <- synteny_df$Start_new[synteny_df$Start_2<6664979]
synteny_df$End_1[synteny_df$Start_2<6664979] <- synteny_df$End_new[synteny_df$Start_2<6664979]

synteny_df <- select(synteny_df, -Start_new, -End_new)


ideogram(karyotype = karyotype_df, synteny = synteny_df,output="chromosome_015_inversion_flipped2.svg")
convertSVG("chromosome_015_inversion_flipped2.svg", file="chromosome_015_inversion_flipped2", device = "png")




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

#karyotype_df$Start[karyotype_df$Chr=="priptg000015l"] <- 17909547
#karyotype_df$End[karyotype_df$Chr=="priptg000015l"] <- 1

#paf_df$qname <- as.numeric(str_remove_all(paf_df$qname,"scaffold_|_RagTag"))
#paf_df$tname <- as.numeric(str_remove_all(paf_df$tname,"scaffold_|_RagTag"))

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

synteny_df$fill[synteny_df$Start_2<12834860] <- "0000FF"
synteny_df$Start_1[synteny_df$Species_1==1] <- 25293088 - synteny_df$Start_1[synteny_df$Species_1==1] 
synteny_df$End_1[synteny_df$Species_1==1] <- 25293088 - synteny_df$End_1[synteny_df$Species_1==1] 

synteny_df$Start_1[synteny_df$Species_1==2] <- 30612747 - synteny_df$Start_1[synteny_df$Species_1==2] 
synteny_df$End_1[synteny_df$Species_1==2] <- 30612747 - synteny_df$End_1[synteny_df$Species_1==2] 


synteny_df$Start_new <- synteny_df$End_1
synteny_df$End_new <- synteny_df$Start_1

synteny_df$Start_1 <- synteny_df$Start_new
synteny_df$End_1 <- synteny_df$End_new


synteny_df <- select(synteny_df, -Start_new, -End_new)





ideogram(karyotype = karyotype_df, synteny = synteny_df,output="chromosome_045_fusion.svg")
convertSVG("chromosome_045_fusion.svg", file="chromosome_045_fusion", device = "png")




###################################################################################
# look at potential fusion!!!!!!

paf_df <- read_paf("scaffold_8_RagTag_8_12_15.paf")

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

synteny_df$fill[synteny_df$Start_2<12834860] <- "0000FF"


ideogram(karyotype = karyotype_df, synteny = synteny_df,output="scaffold_8_RagTag_fusion.svg")
convertSVG("scaffold_8_RagTag_fusion.svg", file="scaffold_8_RagTag_fusion", device = "png")


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

#karyotype_df$Start[karyotype_df$Chr=="priptg000015l"] <- 17909547
#karyotype_df$End[karyotype_df$Chr=="priptg000015l"] <- 1

#paf_df$qname <- as.numeric(str_remove_all(paf_df$qname,"scaffold_|_RagTag"))
#paf_df$tname <- as.numeric(str_remove_all(paf_df$tname,"scaffold_|_RagTag"))

paf_df$qname <- as.numeric(as.factor(paf_df$qname))
paf_df$tname <- as.numeric(as.factor(paf_df$tname))


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


synteny_df$Start_new <- synteny_df$End_1
synteny_df$End_new <- synteny_df$Start_1

synteny_df$Start_1 <- synteny_df$Start_new
synteny_df$End_1 <- synteny_df$End_new


synteny_df <- select(synteny_df, -Start_new, -End_new)



ideogram(karyotype = karyotype_df, synteny = synteny_df,output="prtg00003l_scaffold_13_15_fusion.svg")
convertSVG("prtg00003l_scaffold_13_15_fusion.svg", file="prtg00003l_scaffold_13_15_fusion", device = "png")



















































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





