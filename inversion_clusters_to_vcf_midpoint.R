# inversion genotypes to vcf dummy file of pseudo SNPs
# S. M. Orzechowski, using o3-mini-high for coding guidance
# February 2025


library(dplyr)
library(stringr)
library(ggplot2)

setwd("~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/local_pcangsd/tables")
pops <- read.csv('~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/all_samples_pops_lat_long.csv',header=T)
names <- read.table('~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/LEA/vcf_sample_names_65ind_autos.txt',header=F)

#####################################################################################
table <- read.table('scaffold_4_pop_cluster1_labels.tsv',header=T)
combined <- left_join(pops,table)

# subset to 65 individuals used in the LEA analyses for autosomes
individuals <- str_remove_all(names$V1,"/n/netscratch/edwards_lab/Lab/smorzechowski/meliphagid/analysis/2024-11-03/06-bwa/indel_final/|_bwa_merge_dedup_clip_sort_fixmate_realigned_final.bam|_bwa_dedup_clip_sort_fixmate_realigned_final.bam")

combined_filt <- combined[combined$Sample%in%individuals,]
dummy <- cbind(names,combined_filt)

# Define file name
output_file <- "inversion_dummy_midpoint.vcf"


# Cluster Centroids
#0.130611:homozygote
#-0.127780:homozygote
#-0.003876:heterozygote


# Define genotype mapping as a named vector
genotype_map <- c("0" = "0/0", "2" = "0/1", "1" = "1/1")

#####################################################
data <- dummy
data$Cluster[data$Cluster=="0"] <-"0/0"
data$Cluster[data$Cluster=="2"] <-"0/1"
data$Cluster[data$Cluster=="1"] <-"1/1"

data$sex <- str_extract(data$Sample,"[MFU]$")
data$sex[data$sex=="U"]<-"M"
data$sex[data$Sample=="Nom_093_F"] <- "M"
data$sex[data$Sample=="Nom_094_M"] <- "F"
data$sex <- factor(data$sex)

ggplot(data,aes(Longitude,Cluster,fill=sex,color=sex))+
  geom_point()
#####################################################

# Extract sample names and inversion values
samples <- dummy$V1

# Map inversion values to genotypes (if a value is missing, set to "./.")
genotypes <- genotype_map[as.character(dummy$Cluster)]
genotypes[is.na(genotypes)] <- "./."

# Prepare VCF header lines
header_lines <- c(
  "##fileformat=VCFv4.2",
  "##source=InversionDummyConversion",
  "##INFO=<ID=INV,Number=1,Type=Integer,Description=\"Dummy inversion indicator\">",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
)

# Construct the column header line with sample names
column_header <- paste0(
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
  paste(samples, collapse = "\t")
)

# Create a single variant record using arbitrary values:
#   CHROM = "chrInv" (a pseudo-chromosome for the inversion marker)
#   POS   = 1000     (an arbitrary position)
#   ID    = "inv_dummy"
#   REF   = "A"
#   ALT   = "T"
#   QUAL  = "."
#   FILTER= "PASS"
#   INFO  = "INV=1" (optional, can be used to denote the inversion)
#   FORMAT= "GT"

# Midpoint for inversion 1 on Chr 4: 100.245 Mb

Chr4.1_record_line <- paste(
  "scaffold_4_RagTag",      # CHROM
  "100245000",        # POS
  "inv4.1_dummy",   # ID
  "A",           # REF
  "T",           # ALT
  ".",           # QUAL
  "PASS",        # FILTER
  "INV=1",       # INFO
  "GT",          # FORMAT
  paste(genotypes, collapse = "\t"), # Genotype values for each sample
  sep = "\t"
)


#####################################################################################
table <- read.table('scaffold_4_pop_cluster2_labels.tsv',header=T)
combined <- left_join(pops,table)

# subset to 65 individuals used in the LEA analyses for autosomes
individuals <- str_remove_all(names$V1,"/n/netscratch/edwards_lab/Lab/smorzechowski/meliphagid/analysis/2024-11-03/06-bwa/indel_final/|_bwa_merge_dedup_clip_sort_fixmate_realigned_final.bam|_bwa_dedup_clip_sort_fixmate_realigned_final.bam")

combined_filt <- combined[combined$Sample%in%individuals,]
dummy <- cbind(names,combined_filt)

# Define file name
output_file <- "inversion_dummy_midpoint.vcf"

# Cluster Centroids
#0.119990:homozygote
#-0.166926:homozygote
#-0.035526:heterozygote

# Define genotype mapping as a named vector
genotype_map <- c("0" = "0/0", "2" = "0/1", "1" = "1/1")

#####################################################
data <- dummy
data$Cluster[data$Cluster=="0"] <-"0/0"
data$Cluster[data$Cluster=="2"] <-"0/1"
data$Cluster[data$Cluster=="1"] <-"1/1"

data$sex <- str_extract(data$Sample,"[MFU]$")
data$sex[data$sex=="U"]<-"M"
data$sex[data$Sample=="Nom_093_F"] <- "M"
data$sex[data$Sample=="Nom_094_M"] <- "F"
data$sex <- factor(data$sex)

ggplot(data,aes(Longitude,Cluster,fill=sex,color=sex))+
  geom_point()
#####################################################


# Extract sample names and inversion values
samples <- dummy$V1

# Map inversion values to genotypes (if a value is missing, set to "./.")
genotypes <- genotype_map[as.character(dummy$Cluster)]
genotypes[is.na(genotypes)] <- "./."

# Prepare VCF header lines
header_lines <- c(
  "##fileformat=VCFv4.2",
  "##source=InversionDummyConversion",
  "##INFO=<ID=INV,Number=1,Type=Integer,Description=\"Dummy inversion indicator\">",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
)

# Construct the column header line with sample names
column_header <- paste0(
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
  paste(samples, collapse = "\t")
)

# Create a single variant record using arbitrary values:
#   CHROM = "chrInv" (a pseudo-chromosome for the inversion marker)
#   POS   = 1000     (an arbitrary position)
#   ID    = "inv_dummy"
#   REF   = "A"
#   ALT   = "T"
#   QUAL  = "."
#   FILTER= "PASS"
#   INFO  = "INV=1" (optional, can be used to denote the inversion)
#   FORMAT= "GT"

# Midpoint for inversion 2 on Chr 4: 114.129 Mb

Chr4.2_record_line <- paste(
  "scaffold_4_RagTag",      # CHROM
  "114129000",        # POS
  "inv4.2_dummy",   # ID
  "A",           # REF
  "T",           # ALT
  ".",           # QUAL
  "PASS",        # FILTER
  "INV=1",       # INFO
  "GT",          # FORMAT
  paste(genotypes, collapse = "\t"), # Genotype values for each sample
  sep = "\t"
)


#####################################################################################
table <- read.table('scaffold_6_pop_cluster_labels.tsv',header=T)
combined <- left_join(pops,table)

# subset to 65 individuals used in the LEA analyses for autosomes
individuals <- str_remove_all(names$V1,"/n/netscratch/edwards_lab/Lab/smorzechowski/meliphagid/analysis/2024-11-03/06-bwa/indel_final/|_bwa_merge_dedup_clip_sort_fixmate_realigned_final.bam|_bwa_dedup_clip_sort_fixmate_realigned_final.bam")

combined_filt <- combined[combined$Sample%in%individuals,]
dummy <- cbind(names,combined_filt)

# Define file name
output_file <- "inversion_dummy_midpoint.vcf"

# Cluster Centroids
#0.138682:homozygote
#-0.119128:homozygote
#0.002137:heterozygote


# Define genotype mapping as a named vector
genotype_map <- c("0" = "0/0", "2" = "0/1", "1" = "1/1")


#####################################################
data <- dummy
data$Cluster[data$Cluster=="0"] <-"0/0"
data$Cluster[data$Cluster=="2"] <-"0/1"
data$Cluster[data$Cluster=="1"] <-"1/1"

data$sex <- str_extract(data$Sample,"[MFU]$")
data$sex[data$sex=="U"]<-"M"
data$sex[data$Sample=="Nom_093_F"] <- "M"
data$sex[data$Sample=="Nom_094_M"] <- "F"
data$sex <- factor(data$sex)

ggplot(data,aes(Longitude,Cluster,fill=sex,color=sex))+
  geom_point()
#####################################################

# Extract sample names and inversion values
samples <- dummy$V1

# Map inversion values to genotypes (if a value is missing, set to "./.")
genotypes <- genotype_map[as.character(dummy$Cluster)]
genotypes[is.na(genotypes)] <- "./."

# Prepare VCF header lines
header_lines <- c(
  "##fileformat=VCFv4.2",
  "##source=InversionDummyConversion",
  "##INFO=<ID=INV,Number=1,Type=Integer,Description=\"Dummy inversion indicator\">",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
)

# Construct the column header line with sample names
column_header <- paste0(
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
  paste(samples, collapse = "\t")
)

# Create a single variant record using arbitrary values:
#   CHROM = "chrInv" (a pseudo-chromosome for the inversion marker)
#   POS   = 1000     (an arbitrary position)
#   ID    = "inv_dummy"
#   REF   = "A"
#   ALT   = "T"
#   QUAL  = "."
#   FILTER= "PASS"
#   INFO  = "INV=1" (optional, can be used to denote the inversion)
#   FORMAT= "GT"

# Midpoint for inversion on Chr 6: 24.395 Mb

Chr6_record_line <- paste(
  "scaffold_6_RagTag",      # CHROM
  "24395000",        # POS
  "inv6_dummy",   # ID
  "A",           # REF
  "T",           # ALT
  ".",           # QUAL
  "PASS",        # FILTER
  "INV=1",       # INFO
  "GT",          # FORMAT
  paste(genotypes, collapse = "\t"), # Genotype values for each sample
  sep = "\t"
)

#####################################################################################
table <- read.table('scaffold_8_pop_cluster1_labels.tsv',header=T)
combined <- left_join(pops,table)

# subset to 65 individuals used in the LEA analyses for autosomes
individuals <- str_remove_all(names$V1,"/n/netscratch/edwards_lab/Lab/smorzechowski/meliphagid/analysis/2024-11-03/06-bwa/indel_final/|_bwa_merge_dedup_clip_sort_fixmate_realigned_final.bam|_bwa_dedup_clip_sort_fixmate_realigned_final.bam")

combined_filt <- combined[combined$Sample%in%individuals,]
dummy <- cbind(names,combined_filt)

# Define file name
output_file <- "inversion_dummy_midpoint.vcf"

# Cluster Centroids
#-0.092788:homozygote
#0.054039:heterozygote
#0.191118:homozygote


# Define genotype mapping as a named vector
genotype_map <- c("0" = "0/0", "1" = "0/1", "2" = "1/1")

# Extract sample names and inversion values
samples <- dummy$V1

# Map inversion values to genotypes (if a value is missing, set to "./.")
genotypes <- genotype_map[as.character(dummy$Cluster)]
genotypes[is.na(genotypes)] <- "./."

# Prepare VCF header lines
header_lines <- c(
  "##fileformat=VCFv4.2",
  "##source=InversionDummyConversion",
  "##INFO=<ID=INV,Number=1,Type=Integer,Description=\"Dummy inversion indicator\">",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
)

# Construct the column header line with sample names
column_header <- paste0(
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
  paste(samples, collapse = "\t")
)

# Create a single variant record using arbitrary values:
#   CHROM = "chrInv" (a pseudo-chromosome for the inversion marker)
#   POS   = 1000     (an arbitrary position)
#   ID    = "inv_dummy"
#   REF   = "A"
#   ALT   = "T"
#   QUAL  = "."
#   FILTER= "PASS"
#   INFO  = "INV=1" (optional, can be used to denote the inversion)
#   FORMAT= "GT"

# Midpoint for inversion 1 on Chr 8: ~46.92 Mb
Chr8.1_record_line <- paste(
  "scaffold_8_RagTag",      # CHROM
  "46920000",        # POS
  "inv8.1_dummy",   # ID
  "A",           # REF
  "T",           # ALT
  ".",           # QUAL
  "PASS",        # FILTER
  "INV=1",       # INFO
  "GT",          # FORMAT
  paste(genotypes, collapse = "\t"), # Genotype values for each sample
  sep = "\t"
)

#####################################################################################
table <- read.table('scaffold_8_pop_cluster2_labels.tsv',header=T)
combined <- left_join(pops,table)

# subset to 65 individuals used in the LEA analyses for autosomes
individuals <- str_remove_all(names$V1,"/n/netscratch/edwards_lab/Lab/smorzechowski/meliphagid/analysis/2024-11-03/06-bwa/indel_final/|_bwa_merge_dedup_clip_sort_fixmate_realigned_final.bam|_bwa_dedup_clip_sort_fixmate_realigned_final.bam")

combined_filt <- combined[combined$Sample%in%individuals,]
dummy <- cbind(names,combined_filt)

# Define file name
output_file <- "inversion_dummy_midpoint.vcf"

# Cluster Centroids
#-0.119055:homozygote
#0.164138:homozygote
#0.037786:heterozygote


# Define genotype mapping as a named vector
genotype_map <- c("0" = "0/0", "2" = "0/1", "1" = "1/1")

# Extract sample names and inversion values
samples <- dummy$V1

# Map inversion values to genotypes (if a value is missing, set to "./.")
genotypes <- genotype_map[as.character(dummy$Cluster)]
genotypes[is.na(genotypes)] <- "./."

# Prepare VCF header lines
header_lines <- c(
  "##fileformat=VCFv4.2",
  "##source=InversionDummyConversion",
  "##INFO=<ID=INV,Number=1,Type=Integer,Description=\"Dummy inversion indicator\">",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
)

# Construct the column header line with sample names
column_header <- paste0(
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
  paste(samples, collapse = "\t")
)

# Create a single variant record using arbitrary values:
#   CHROM = "chrInv" (a pseudo-chromosome for the inversion marker)
#   POS   = 1000     (an arbitrary position)
#   ID    = "inv_dummy"
#   REF   = "A"
#   ALT   = "T"
#   QUAL  = "."
#   FILTER= "PASS"
#   INFO  = "INV=1" (optional, can be used to denote the inversion)
#   FORMAT= "GT"

# Midpoint for inversion 2 on Chr 8: ~ 5.36 Mb

Chr8.2_record_line <- paste(
  "scaffold_8_RagTag",      # CHROM
  "5360000",        # POS
  "inv8.2_dummy",   # ID
  "A",           # REF
  "T",           # ALT
  ".",           # QUAL
  "PASS",        # FILTER
  "INV=1",       # INFO
  "GT",          # FORMAT
  paste(genotypes, collapse = "\t"), # Genotype values for each sample
  sep = "\t"
)


#####################################################################################
table <- read.table('scaffold_9_pop_cluster_labels.tsv',header=T)
combined <- left_join(pops,table)

# subset to 65 individuals used in the LEA analyses for autosomes
individuals <- str_remove_all(names$V1,"/n/netscratch/edwards_lab/Lab/smorzechowski/meliphagid/analysis/2024-11-03/06-bwa/indel_final/|_bwa_merge_dedup_clip_sort_fixmate_realigned_final.bam|_bwa_dedup_clip_sort_fixmate_realigned_final.bam")

combined_filt <- combined[combined$Sample%in%individuals,]
dummy <- cbind(names,combined_filt)

# Define file name
output_file <- "inversion_dummy_midpoint.vcf"

# Cluster Centroids
#-0.029322:heterozygote
#0.097375:homozygote
#-0.168284:homozygote


# Define genotype mapping as a named vector
genotype_map <- c("1" = "0/0", "0" = "0/1", "2" = "1/1")

# Extract sample names and inversion values
samples <- dummy$V1

# Map inversion values to genotypes (if a value is missing, set to "./.")
genotypes <- genotype_map[as.character(dummy$Cluster)]
genotypes[is.na(genotypes)] <- "./."

# Prepare VCF header lines
header_lines <- c(
  "##fileformat=VCFv4.2",
  "##source=InversionDummyConversion",
  "##INFO=<ID=INV,Number=1,Type=Integer,Description=\"Dummy inversion indicator\">",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
)

# Construct the column header line with sample names
column_header <- paste0(
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
  paste(samples, collapse = "\t")
)

# Create a single variant record using arbitrary values:
#   CHROM = "chrInv" (a pseudo-chromosome for the inversion marker)
#   POS   = 1000     (an arbitrary position)
#   ID    = "inv_dummy"
#   REF   = "A"
#   ALT   = "T"
#   QUAL  = "."
#   FILTER= "PASS"
#   INFO  = "INV=1" (optional, can be used to denote the inversion)
#   FORMAT= "GT"

# Midpoint for inversion on Chr 9: ~ 9.98 Mb 

Chr9_record_line <- paste(
  "scaffold_9_RagTag",      # CHROM
  "9980000",        # POS
  "inv9_dummy",   # ID
  "A",           # REF
  "T",           # ALT
  ".",           # QUAL
  "PASS",        # FILTER
  "INV=1",       # INFO
  "GT",          # FORMAT
  paste(genotypes, collapse = "\t"), # Genotype values for each sample
  sep = "\t"
)

#####################################################################################
table <- read.table('scaffold_10_pop_cluster_labels.tsv',header=T)
combined <- left_join(pops,table)

# subset to 65 individuals used in the LEA analyses for autosomes
individuals <- str_remove_all(names$V1,"/n/netscratch/edwards_lab/Lab/smorzechowski/meliphagid/analysis/2024-11-03/06-bwa/indel_final/|_bwa_merge_dedup_clip_sort_fixmate_realigned_final.bam|_bwa_dedup_clip_sort_fixmate_realigned_final.bam")

combined_filt <- combined[combined$Sample%in%individuals,]
dummy <- cbind(names,combined_filt)

# Define file name
output_file <- "inversion_dummy_midpoint.vcf"

# Cluster Centroids
#-0.123191:homozygote
#0.136582:homozygote
#0.010261:heterozygote


# Define genotype mapping as a named vector
genotype_map <- c("0" = "0/0", "2" = "0/1", "1" = "1/1")

# Extract sample names and inversion values
samples <- dummy$V1

# Map inversion values to genotypes (if a value is missing, set to "./.")
genotypes <- genotype_map[as.character(dummy$Cluster)]
genotypes[is.na(genotypes)] <- "./."

# Prepare VCF header lines
header_lines <- c(
  "##fileformat=VCFv4.2",
  "##source=InversionDummyConversion",
  "##INFO=<ID=INV,Number=1,Type=Integer,Description=\"Dummy inversion indicator\">",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
)

# Construct the column header line with sample names
column_header <- paste0(
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
  paste(samples, collapse = "\t")
)

# Create a single variant record using arbitrary values:
#   CHROM = "chrInv" (a pseudo-chromosome for the inversion marker)
#   POS   = 1000     (an arbitrary position)
#   ID    = "inv_dummy"
#   REF   = "A"
#   ALT   = "T"
#   QUAL  = "."
#   FILTER= "PASS"
#   INFO  = "INV=1" (optional, can be used to denote the inversion)
#   FORMAT= "GT"

# Midpoint for inversion on Chr 10: 16.005 Mb
Chr10_record_line <- paste(
  "scaffold_10_RagTag",      # CHROM
  "16005000",        # POS
  "inv10_dummy",   # ID
  "A",           # REF
  "T",           # ALT
  ".",           # QUAL
  "PASS",        # FILTER
  "INV=1",       # INFO
  "GT",          # FORMAT
  paste(genotypes, collapse = "\t"), # Genotype values for each sample
  sep = "\t"
)


#####################################################################################
table <- read.table('scaffold_13_pop_cluster_labels.tsv',header=T)
combined <- left_join(pops,table)

# subset to 65 individuals used in the LEA analyses for autosomes
individuals <- str_remove_all(names$V1,"/n/netscratch/edwards_lab/Lab/smorzechowski/meliphagid/analysis/2024-11-03/06-bwa/indel_final/|_bwa_merge_dedup_clip_sort_fixmate_realigned_final.bam|_bwa_dedup_clip_sort_fixmate_realigned_final.bam")

combined_filt <- combined[combined$Sample%in%individuals,]
dummy <- cbind(names,combined_filt)

# Define file name
output_file <- "inversion_dummy_midpoint.vcf"

# Cluster Centroids
#0.102898:homozygote
#-0.170139:homozygote
#-0.032086:heterozygote


# Define genotype mapping as a named vector
genotype_map <- c("0" = "0/0", "2" = "0/1", "1" = "1/1")

# Extract sample names and inversion values
samples <- dummy$V1

# Map inversion values to genotypes (if a value is missing, set to "./.")
genotypes <- genotype_map[as.character(dummy$Cluster)]
genotypes[is.na(genotypes)] <- "./."

# Prepare VCF header lines
header_lines <- c(
  "##fileformat=VCFv4.2",
  "##source=InversionDummyConversion",
  "##INFO=<ID=INV,Number=1,Type=Integer,Description=\"Dummy inversion indicator\">",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
)

# Construct the column header line with sample names
column_header <- paste0(
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
  paste(samples, collapse = "\t")
)

# Create a single variant record using arbitrary values:
#   CHROM = "chrInv" (a pseudo-chromosome for the inversion marker)
#   POS   = 1000     (an arbitrary position)
#   ID    = "inv_dummy"
#   REF   = "A"
#   ALT   = "T"
#   QUAL  = "."
#   FILTER= "PASS"
#   INFO  = "INV=1" (optional, can be used to denote the inversion)
#   FORMAT= "GT"

# Midpoint for inversion on Chr 13: ~ 20.4 Mb
Chr13_record_line <- paste(
  "scaffold_13_RagTag",      # CHROM
  "20400000",        # POS
  "inv13_dummy",   # ID
  "A",           # REF
  "T",           # ALT
  ".",           # QUAL
  "PASS",        # FILTER
  "INV=1",       # INFO
  "GT",          # FORMAT
  paste(genotypes, collapse = "\t"), # Genotype values for each sample
  sep = "\t"
)


#####################################################################################
table <- read.table('scaffold_14_pop_cluster_labels.tsv',header=T)
combined <- left_join(pops,table)

# subset to 65 individuals used in the LEA analyses for autosomes
individuals <- str_remove_all(names$V1,"/n/netscratch/edwards_lab/Lab/smorzechowski/meliphagid/analysis/2024-11-03/06-bwa/indel_final/|_bwa_merge_dedup_clip_sort_fixmate_realigned_final.bam|_bwa_dedup_clip_sort_fixmate_realigned_final.bam")

combined_filt <- combined[combined$Sample%in%individuals,]
dummy <- cbind(names,combined_filt)

# Define file name
output_file <- "inversion_dummy_midpoint.vcf"

# Cluster Centroids
#0.099481:homozygote
#-0.170824:homozygote
#-0.037793:heterozygote


# Define genotype mapping as a named vector
genotype_map <- c("0" = "0/0", "2" = "0/1", "1" = "1/1")

# Extract sample names and inversion values
samples <- dummy$V1

# Map inversion values to genotypes (if a value is missing, set to "./.")
genotypes <- genotype_map[as.character(dummy$Cluster)]
genotypes[is.na(genotypes)] <- "./."

# Prepare VCF header lines
header_lines <- c(
  "##fileformat=VCFv4.2",
  "##source=InversionDummyConversion",
  "##INFO=<ID=INV,Number=1,Type=Integer,Description=\"Dummy inversion indicator\">",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
)

# Construct the column header line with sample names
column_header <- paste0(
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
  paste(samples, collapse = "\t")
)

# Create a single variant record using arbitrary values:
#   CHROM = "chrInv" (a pseudo-chromosome for the inversion marker)
#   POS   = 1000     (an arbitrary position)
#   ID    = "inv_dummy"
#   REF   = "A"
#   ALT   = "T"
#   QUAL  = "."
#   FILTER= "PASS"
#   INFO  = "INV=1" (optional, can be used to denote the inversion)
#   FORMAT= "GT"

# Midpoint for inversion on Chr 14: 5.085 Mb
Chr14_record_line <- paste(
  "scaffold_14_RagTag",      # CHROM
  "5085000",        # POS
  "inv14_dummy",   # ID
  "A",           # REF
  "T",           # ALT
  ".",           # QUAL
  "PASS",        # FILTER
  "INV=1",       # INFO
  "GT",          # FORMAT
  paste(genotypes, collapse = "\t"), # Genotype values for each sample
  sep = "\t"
)

#####################################################################################
table <- read.table('scaffold_15pop_cluster_labels.tsv',header=T)
combined <- left_join(pops,table)

# subset to 65 individuals used in the LEA analyses for autosomes
individuals <- str_remove_all(names$V1,"/n/netscratch/edwards_lab/Lab/smorzechowski/meliphagid/analysis/2024-11-03/06-bwa/indel_final/|_bwa_merge_dedup_clip_sort_fixmate_realigned_final.bam|_bwa_dedup_clip_sort_fixmate_realigned_final.bam")

combined_filt <- combined[combined$Sample%in%individuals,]
dummy <- cbind(names,combined_filt)

# Define file name
output_file <- "inversion_dummy_midpoint.vcf"

# Cluster Centroids
#0.104248:homozygote
#-0.167038:homozygote
#-0.037992:heterozygote


# Define genotype mapping as a named vector
genotype_map <- c("0" = "0/0", "2" = "0/1", "1" = "1/1")

# Extract sample names and inversion values
samples <- dummy$V1

# Map inversion values to genotypes (if a value is missing, set to "./.")
genotypes <- genotype_map[as.character(dummy$Cluster)]
genotypes[is.na(genotypes)] <- "./."

# Prepare VCF header lines
header_lines <- c(
  "##fileformat=VCFv4.2",
  "##source=InversionDummyConversion",
  "##INFO=<ID=INV,Number=1,Type=Integer,Description=\"Dummy inversion indicator\">",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
)

# Construct the column header line with sample names
column_header <- paste0(
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
  paste(samples, collapse = "\t")
)

# Create a single variant record using arbitrary values:
#   CHROM = "chrInv" (a pseudo-chromosome for the inversion marker)
#   POS   = 1000     (an arbitrary position)
#   ID    = "inv_dummy"
#   REF   = "A"
#   ALT   = "T"
#   QUAL  = "."
#   FILTER= "PASS"
#   INFO  = "INV=1" (optional, can be used to denote the inversion)
#   FORMAT= "GT"

# Midpoint for inversion on Chr 15: 4.385 Mb
Chr15_record_line <- paste(
  "scaffold_15_RagTag",      # CHROM
  "4385000",        # POS
  "inv15_dummy",   # ID
  "A",           # REF
  "T",           # ALT
  ".",           # QUAL
  "PASS",        # FILTER
  "INV=1",       # INFO
  "GT",          # FORMAT
  paste(genotypes, collapse = "\t"), # Genotype values for each sample
  sep = "\t"
)


#####################################################################################
table <- read.table('scaffold_17_pop_cluster_labels.tsv',header=T)
combined <- left_join(pops,table)

# subset to 65 individuals used in the LEA analyses for autosomes
individuals <- str_remove_all(names$V1,"/n/netscratch/edwards_lab/Lab/smorzechowski/meliphagid/analysis/2024-11-03/06-bwa/indel_final/|_bwa_merge_dedup_clip_sort_fixmate_realigned_final.bam|_bwa_dedup_clip_sort_fixmate_realigned_final.bam")

combined_filt <- combined[combined$Sample%in%individuals,]
dummy <- cbind(names,combined_filt)

# Define file name
output_file <- "inversion_dummy_midpoint.vcf"

# Cluster Centroids
#0.189696:homozygotes
#-0.094570:homozygotes
#0.063400:heterozygotes

# Define genotype mapping as a named vector
genotype_map <- c("0" = "0/0", "2" = "0/1", "1" = "1/1")

# Extract sample names and inversion values
samples <- dummy$V1

# Map inversion values to genotypes (if a value is missing, set to "./.")
genotypes <- genotype_map[as.character(dummy$Cluster)]
genotypes[is.na(genotypes)] <- "./."

# Prepare VCF header lines
header_lines <- c(
  "##fileformat=VCFv4.2",
  "##source=InversionDummyConversion",
  "##INFO=<ID=INV,Number=1,Type=Integer,Description=\"Dummy inversion indicator\">",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
)

# Construct the column header line with sample names
column_header <- paste0(
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
  paste(samples, collapse = "\t")
)

# Create a single variant record using arbitrary values:
#   CHROM = "chrInv" (a pseudo-chromosome for the inversion marker)
#   POS   = 1000     (an arbitrary position)
#   ID    = "inv_dummy"
#   REF   = "A"
#   ALT   = "T"
#   QUAL  = "."
#   FILTER= "PASS"
#   INFO  = "INV=1" (optional, can be used to denote the inversion)
#   FORMAT= "GT"

# Midpoint for inversion on Chr 17: 16.27 Mb
Chr17_record_line <- paste(
  "scaffold_17_RagTag",      # CHROM
  "16270000",        # POS
  "inv17_dummy",   # ID
  "A",           # REF
  "T",           # ALT
  ".",           # QUAL
  "PASS",        # FILTER
  "INV=1",       # INFO
  "GT",          # FORMAT
  paste(genotypes, collapse = "\t"), # Genotype values for each sample
  sep = "\t"
)



#####################################################################################
table <- read.table('scaffold_18_pop_cluster_labels.tsv',header=T)
combined <- left_join(pops,table)

# subset to 65 individuals used in the LEA analyses for autosomes
individuals <- str_remove_all(names$V1,"/n/netscratch/edwards_lab/Lab/smorzechowski/meliphagid/analysis/2024-11-03/06-bwa/indel_final/|_bwa_merge_dedup_clip_sort_fixmate_realigned_final.bam|_bwa_dedup_clip_sort_fixmate_realigned_final.bam")

combined_filt <- combined[combined$Sample%in%individuals,]
dummy <- cbind(names,combined_filt)

# Define file name
output_file <- "inversion_dummy_midpoint.vcf"

# Cluster Centroids
#0.132701:homozygotes
#-0.132051:homozygotes
#0.001952:heterozygotes

# Define genotype mapping as a named vector
genotype_map <- c("0" = "0/0", "2" = "0/1", "1" = "1/1")

# Extract sample names and inversion values
samples <- dummy$V1

# Map inversion values to genotypes (if a value is missing, set to "./.")
genotypes <- genotype_map[as.character(dummy$Cluster)]
genotypes[is.na(genotypes)] <- "./."

# Prepare VCF header lines
header_lines <- c(
  "##fileformat=VCFv4.2",
  "##source=InversionDummyConversion",
  "##INFO=<ID=INV,Number=1,Type=Integer,Description=\"Dummy inversion indicator\">",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
)

# Construct the column header line with sample names
column_header <- paste0(
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
  paste(samples, collapse = "\t")
)

# Create a single variant record using arbitrary values:
#   CHROM = "chrInv" (a pseudo-chromosome for the inversion marker)
#   POS   = 1000     (an arbitrary position)
#   ID    = "inv_dummy"
#   REF   = "A"
#   ALT   = "T"
#   QUAL  = "."
#   FILTER= "PASS"
#   INFO  = "INV=1" (optional, can be used to denote the inversion)
#   FORMAT= "GT"

# Midpoint for inversion on Chr 18: 6.17 Mb

Chr18_record_line <- paste(
  "scaffold_18_RagTag",      # CHROM
  "6170000",        # POS
  "inv18_dummy",   # ID
  "A",           # REF
  "T",           # ALT
  ".",           # QUAL
  "PASS",        # FILTER
  "INV=1",       # INFO
  "GT",          # FORMAT
  paste(genotypes, collapse = "\t"), # Genotype values for each sample
  sep = "\t"
)

#####################################################################################
table <- read.table('pri_ptg000045l_pop_cluster_labels.tsv',header=T)
combined <- left_join(pops,table)

# subset to 65 individuals used in the LEA analyses for autosomes
individuals <- str_remove_all(names$V1,"/n/netscratch/edwards_lab/Lab/smorzechowski/meliphagid/analysis/2024-11-03/06-bwa/indel_final/|_bwa_merge_dedup_clip_sort_fixmate_realigned_final.bam|_bwa_dedup_clip_sort_fixmate_realigned_final.bam")

combined_filt <- combined[combined$Sample%in%individuals,]
dummy <- cbind(names,combined_filt)

# Define file name
output_file <- "inversion_dummy_midpoint.vcf"

# Cluster Centroids
#0.120864:homozygotes
#-0.137718:homozygotes
#-0.012890:heterozygotes

# Define genotype mapping as a named vector
genotype_map <- c("0" = "0/0", "2" = "0/1", "1" = "1/1")

# Extract sample names and inversion values
samples <- dummy$V1

# Map inversion values to genotypes (if a value is missing, set to "./.")
genotypes <- genotype_map[as.character(dummy$Cluster)]
genotypes[is.na(genotypes)] <- "./."

# Prepare VCF header lines
header_lines <- c(
  "##fileformat=VCFv4.2",
  "##source=InversionDummyConversion",
  "##INFO=<ID=INV,Number=1,Type=Integer,Description=\"Dummy inversion indicator\">",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
)

# Construct the column header line with sample names
column_header <- paste0(
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
  paste(samples, collapse = "\t")
)

# Create a single variant record using arbitrary values:
#   CHROM = "chrInv" (a pseudo-chromosome for the inversion marker)
#   POS   = 1000     (an arbitrary position)
#   ID    = "inv_dummy"
#   REF   = "A"
#   ALT   = "T"
#   QUAL  = "."
#   FILTER= "PASS"
#   INFO  = "INV=1" (optional, can be used to denote the inversion)
#   FORMAT= "GT"

# Midpoint for inversion on scaffold 45: 6.415 Mb

ptg45_record_line <- paste(
  "pri#ptg000045l",      # CHROM
  "6415000",        # POS
  "inv_ptg45_dummy",   # ID
  "A",           # REF
  "T",           # ALT
  ".",           # QUAL
  "PASS",        # FILTER
  "INV=1",       # INFO
  "GT",          # FORMAT
  paste(genotypes, collapse = "\t"), # Genotype values for each sample
  sep = "\t"
)

#####################################################################################
table <- read.table('pri_ptg000015l_pop_cluster_labels.tsv',header=T)
combined <- left_join(pops,table)

# subset to 65 individuals used in the LEA analyses for autosomes
individuals <- str_remove_all(names$V1,"/n/netscratch/edwards_lab/Lab/smorzechowski/meliphagid/analysis/2024-11-03/06-bwa/indel_final/|_bwa_merge_dedup_clip_sort_fixmate_realigned_final.bam|_bwa_dedup_clip_sort_fixmate_realigned_final.bam")

combined_filt <- combined[combined$Sample%in%individuals,]
dummy <- cbind(names,combined_filt)

# Define file name
output_file <- "inversion_dummy_midpoint.vcf"

# Cluster Centroids
#-0.139030:homozygotes
#0.128042:homozygotes
#-0.003443:heterozygotes

# Define genotype mapping as a named vector
genotype_map <- c("0" = "0/0", "2" = "0/1", "1" = "1/1")

# Extract sample names and inversion values
samples <- dummy$V1

# Map inversion values to genotypes (if a value is missing, set to "./.")
genotypes <- genotype_map[as.character(dummy$Cluster)]
genotypes[is.na(genotypes)] <- "./."

# Prepare VCF header lines
header_lines <- c(
  "##fileformat=VCFv4.2",
  "##source=InversionDummyConversion",
  "##INFO=<ID=INV,Number=1,Type=Integer,Description=\"Dummy inversion indicator\">",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
)

# Construct the column header line with sample names
column_header <- paste0(
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
  paste(samples, collapse = "\t")
)

# Create a single variant record using arbitrary values:
#   CHROM = "chrInv" (a pseudo-chromosome for the inversion marker)
#   POS   = 1000     (an arbitrary position)
#   ID    = "inv_dummy"
#   REF   = "A"
#   ALT   = "T"
#   QUAL  = "."
#   FILTER= "PASS"
#   INFO  = "INV=1" (optional, can be used to denote the inversion)
#   FORMAT= "GT"

# Midpoint for inversion on scaffold 15: 12.385 Mb

ptg15_record_line <- paste(
  "pri#ptg000015l",      # CHROM
  "12385000",        # POS
  "inv_ptg15_dummy",   # ID
  "A",           # REF
  "T",           # ALT
  ".",           # QUAL
  "PASS",        # FILTER
  "INV=1",       # INFO
  "GT",          # FORMAT
  paste(genotypes, collapse = "\t"), # Genotype values for each sample
  sep = "\t"
)

#####################################################################################
# Write the VCF to the output file
con <- file(output_file, "w")
writeLines(header_lines, con)
writeLines(column_header, con)
writeLines(Chr4.1_record_line, con)
writeLines(Chr4.2_record_line, con)
writeLines(Chr6_record_line, con)
writeLines(Chr8.2_record_line, con)
writeLines(Chr8.1_record_line, con)
writeLines(Chr9_record_line, con)
writeLines(Chr10_record_line, con)
writeLines(Chr13_record_line, con)
writeLines(Chr14_record_line, con)
writeLines(Chr15_record_line, con)
writeLines(Chr17_record_line, con)
writeLines(Chr18_record_line, con)
writeLines(ptg45_record_line, con)
writeLines(ptg15_record_line, con)
close(con)

cat("VCF file written to", output_file, "\n")
