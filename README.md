# Scripts associated with a population genomic analysis of White-eared Honeyeaters (*Nesoptilotis leucotis*)

This readme provides an overview of the genomic analyses conducted in this project. I primarily ran the ANGSD pipeline for variant calling because I sequenced 72 individuals at an average depth of 6x coverage, which is in the low-to-mid coverage range better suited for the genotype likelihood framework of ANGSD. Howevever, I include the manual GATK pipeline for hard calling variants as well.  

The software and programs used in this project include:

- [hifiasm](https://github.com/chhylp123/hifiasm)
- [findZX](https://github.com/hsigeman/findZX)
- [RepeatModeler](https://github.com/Dfam-consortium/RepeatModeler/tree/master)
- [RepeatMasker](https://github.com/Dfam-consortium/RepeatMasker)
- [trimgalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
  - FastQC
- [BWA (v.7.17)](https://bio-bwa.sourceforge.net/bwa.shtml)
- [picard](https://broadinstitute.github.io/picard/)
- [samtools (v1.20)](https://www.htslib.org/doc/1.20/samtools.html) 
- [GATK](https://gatk.broadinstitute.org/hc/en-us)
- [ANGSD](https://www.popgen.dk/angsd/index.php/ANGSD)
- [PCAngsd](https://github.com/Rosemeis/pcangsd)
- [ngsRelate](https://github.com/ANGSD/NgsRelate)
- [vcftools](https://vcftools.github.io/index.html)
- [bcftools](https://samtools.github.io/bcftools/bcftools.html)
- [qualimap](http://qualimap.conesalab.org/)
- [R (v4.3.2)](https://www.r-project.org/)


# Contents
- [Order of operations for each pipeline](#order-of-operations-for-each-pipeline)
- [Genome assembly and curation](#genome-assembly-and-curation)
- [Adapter trimming and read mapping](#adapter-trimming-and-read-mapping)
### [ANGSD pipeline](#angsd-pipeline)
- [Flagging duplicates and clipping reads](#flagging-duplicates-and-clipping-reads)
### [GATK pipeline](#gatk-pipeline)


## Order of operations for each pipeline

### GATK
1. trimgalore.jobscript
2. bwa_sort.jobscript
3. merge_dedup_flag.jobscript
4. split genome into smaller intervals with split_fastaindex.py
5. haplotypecaller.jobscript
6. genomicsDB_import.jobscript
7. genotypeGVCF.jobscript
8. Combine VCF files
9. Filter VCF files
10. Run PCA

### ANGSD
1. trimgalore.jobscript
2. bwa_sort.jobscript
3. merge_dedup_clip.jobscript
4. indel_realign.jobscript


## Genome assembly and curation

I generated HiFi PacBio reads on a Revio SMRTCELL from a female White-eared Honeyeater to create a long-read reference genome with both the neo-Z and neo-W sex chromosomes. I used Heng Li's program `hifiasm` to conduct *de novo* assembly of the HiFi reads. I adapted a snakemake pipeline from Danielle Khost of the Harvard Informatics group (https://github.com/harvardinformatics/pacbio_hifi_assembly).

My config.yaml file for the snakemake pipeline is here:

```
##############################
# Variables you need to change
##############################

sample: "nleucotis_yamma_366917_f"  #Name of the sample to act as base name
reads: "/n/holylfs04/LABS/edwards_lab/Lab/smorzechowski/meliphagid/HiFi_data/1_D01/hifi_reads/m84147_240713_010125_s4.hifi_reads.bam" #List of BAM files output from sequencer (include full path). If multiple files, separate by SPACES

#Hifiasm options
use_hiC: "false"  #Indicate whether using HiC data to aid in with assembly (note, NOT the same as HiC scaffolding!)
hiC_read1: " "  #If using HiC, include path to first set of reads in pair
hiC_read2: " "  ##If using HiC, include path to second set of reads in pair

#Optional (ultra-long) Nanopore integration
use_ont: "false" #Indicate whether have nanopore data
ont_reads: " " #If using Nanopore, add path to COMBINED reads as a fastq file (opt. gzipped fastq)

#Variables you DON'T need to change:
#Output naming:
hifi_outdir: "results/hifiasm_output_RAW"
assem_outdir: "results/assembly"

```

Hifiasm is able to produce two haplotype assemblies (hap1 and hap2) as well as a primary haploid assembly that collapses the bubbles in the graph. Since I am not doing pangenome analysis, I can use the primary assembly to do short read alignment for variant calling. However, it is often necessary to manually curate the sex chromosomes, especially in a species like White-eared Honeyeater, which contains neo-sex chromosomes. I found that the primary assembly did not contain the full phased sequence for the neo-Z or neo-W because the primary assembly algorithm incorrectly purged haplotigs that it characterized as redundant but were actually minimally diverged regions of the neo-sex chromosomes. For instance, it purged the contig(s) containing the new PAR because it is not diverged at all between the neo-sex chromosomes, however it also purged regions where recombination suppression has occurred on the neo-W, leading to divergence between the Z and W gametologs. 

To make sure that all W-linked and Z-linked contigs were included in the final assembly, I had to identify the Z and W linked contigs in all assemblies with a combination of sex differences in coverage and synteny analysis.

I ran findZX to conduct a synteny analysis between White-eared Honeyeater and Blue-faced Honeyeater and estimate coverage of the short reads I aligned across the genome from one male and one female. See the R scripts `classify_contigs_cyan_synteny_hap1.R`, `classify_contigs_cyan_cynteny_hap2.R` and `classify_contigs_tgut_synteny_primary.R` in files. 

Here is an example of how I used the findZX output to find contigs in the White-eared Honeyeater genome that were homologous to the neo-Z (`scaffold_2`) of the Blue-faced Honeyeater genome:

```
setwd("~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE assembly/findZX/Nesoptilotis_hap2")

# get list of autosomal contigs to estimate median depth for each sex
# https://stackoverflow.com/questions/67192258/cant-read-txt-in-r-with-hash-tag-in-the-txt-file

lengths <- read.table("coverage/nleucotis_yamma_366917_f.hap2.p_ctg.filter.10000.fasta.fai",header=F,sep="\t",comment.char="")
matches <- read.table("synteny_lastal/E_cyanotis/bestMatch.list",header=F,sep="\t",comment.char = "")


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
```

And here is an example for the neo-W contigs:

```
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
```

Next, I used the FindZX coverage data for the male and female individual I mapped to the White-eared Honeyeater assemblies (hap1, hap2, primary) to estimate the metric `log2FM` following Schield et al. 2022 (Schield DR, Perry BW, Card DC, Pasquesi GIM, Westfall AK, Mackessy SP, Castoe TA. 2022. The Rattlesnake W Chromosome: A GC-Rich Retroelement Refugium with Retained Gene Function Across Ancient Evolutionary Strata. Genome Biol. Evol. 14:evac116.)

Here's an example:

```
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

# Subset coverage based on window hits to neo-Z
coverage_neoZ <- coverage[coverage$V1%in%matches_neoZ_final$V1,]
coverage_neoZ$contig <- str_remove_all(coverage_neoZ$V1,"nleucotis_yamma_366917_f.hap2#pri#")


ggplot(coverage_neoZ,aes(contig,log2FM))+
  geom_boxplot()

```

Next, I wanted to softmask the repetitive content of the genome to speed up alignment of short reads. Repetitive regions make alignment very difficult and aligners can get stuck in computational black holes trying to align very repetitive short reads to repetitive genome regions that are identical in hundreds or thousands of places across the genome. Variant calling from short reads in repetitive regions is basically hopeless (although I should note that pangenomes allow much better characterization of structural and sequence variants in repetitive regions), thus it is best to focus on the non-repetitive portion of the genome where we can be more confident about the variants we call.


I created a *de novo* repeat library of the White-eared Honeyeater with RepeatModeler:

```
RepeatModeler='/n/home09/smorzechowski/bin/RepeatModeler-2.0.2a' 

# $RepeatModeler/BuildDatabase -name Nleucotis_hifi_v1.0 -engine ncbi $assembly 
$RepeatModeler/RepeatModeler -pa 12 -engine ncbi -database Nleucotis_hifi_v1.0 2>&1 | tee repeatmodeler.log
```

Then I softmasked all repetitive regions using RepeatMasker, supplying the custom RepeatModeler library:

```
RepeatModeler_library=$1
Genome=$2
DIR=$3

export PERL5LIB=/

RepeatMaskerPath='/n/home09/smorzechowski/bin/RepeatMasker'
# softmask = xsmall
mkdir $DIR
$RepeatMaskerPath/RepeatMasker -pa 8 -e ncbi -xsmall -lib $RepeatModeler_library -dir $DIR $Genome
```

Then, I used RagTag to scaffold the contigs into chromosomes based on the Blue-faced Honeyeater, which is the most closely related species of honeyeater with a HiC-scaffolded genome that I produced for a different manuscript:

```
# Scaffold Nleucotis manually-curated softmasked assembly to Blue-faced Honeyeater HiC scaffolded assembly
#query='/n/holylfs04/LABS/edwards_lab/Lab/smorzechowski/meliphagid/ReferenceAssemblies/Nleucotis_hifi_v1.0_headclean_softmasked_fcs_fx.fasta'
#ref='/n/holylfs04/LABS/edwards_lab/Lab/smorzechowski/meliphagid/ReferenceAssemblies/Ecyan_HiC_v3.0.fa'
#outdir='ref_Ecyan_HiC_v3.0_query_Nleu_v1.0'
#species='Nleu'

module load python
source activate ragtag

ref=$1
query=$2
out=$3

ragtag.py scaffold -t 6 $ref $query -r -o $out
```

Then, I also uploaded the genome to NCBI, which runs contamination screening on the genome to see if adapter contamination or sequences from foreign organisms are present. 

I also hardmasked the new PAR on the neo-W to avoid mapping issues which occur when reads map to multiple identical sequences in a genome:

```
# interval information
scaffold_5_RagTag 0 26782956 # added W
scaffold_5_RagTag 26782957 32320175 # ancestral W
scaffold_5_RagTag 32320176 55423182 # added W
scaffold_5_RagTag 55423183 73135152 # ancestral W  
scaffold_5_RagTag 73135153 90514351 # new PAR -- WHICH SHOULD BE HARDMASKED!!!!!
Size of the newPAR based on scaffold_5 = 17,379,198 bp

newPAR.bed
scaffold_5_RagTag 73135153 90514351

mamba activate bedtools

# by default, hard masking is done 
bedtools maskfasta -fi $MEL/ReferenceAssemblies/Nleucotis_hifi_v1.0_hc_sm_fx_scaffolded.fasta -bed newPAR.bed -fo Nleucotis_hifi_v1.0_hc_sm_fx_scaffolded_PAR_masked.fasta


```
The final genome assembly I used in all sequences was: 

`Nleucotis_hifi_v1.0_hc_sm_fx_scaffolded_PAR_masked.fasta`

## Adapter trimming and read mapping

I used trimgalore to trim adapters from the Illumina short reads I generated for 72 individuals of White-eared Honeyeaters. 
See trimgalore.jobscript in files. 

```
module purge
module load python
source activate cutadapt
# this env includes fastqc as well


# $1 = R1 forward reads
# $2 = R2 reverse reads

TRIMGALORE="/n/home09/smorzechowski/bin/TrimGalore-0.6.6/trim_galore"

# SEE OPTIONS $TRIMGALORE --help

# remove retained_unpaired, length=36 for paired
# set quality to 0 to include all base qualities, which ANGSD can handle and filter
# removing bases with quality less than 20 is standard for GATK and other variant callers

$TRIMGALORE --quality 20 --phred33 --fastqc --paired --output_dir trimmed_reads_all_qual --length 36 --stringency 3 -e 0.1 $1 $2
```

To run trimgalore on many individuals, I created a trimgalore.runscript using sed/awk to submit each trimgalore job to SLURM.
See trimgalore.runscript in files. The runscript is a bash script that calls the jobscript and provides filepath to each forward and reverse fastq file.

For example:

```
#!/bin/bash
sbatch trimgalore.jobscript data/P536_SOrze_A01v1_Nom_093_F_S1_L003_R1_001.fastq.gz data/P536_SOrze_A01v1_Nom_093_F_S1_L003_R2_001.fastq.gz
sbatch trimgalore.jobscript data/P536_SOrze_A01v1_Nom_093_F_S1_L004_R1_001.fastq.gz data/P536_SOrze_A01v1_Nom_093_F_S1_L004_R2_001.fastq.gz
```

To map (i.e. align) the trimmed reads to the reference genome, I used BWA to create `.bam` alignment files. 

First I define input parameters
```
# Input parameters
RAW_DIR='/n/holyscratch01/edwards_lab/smorzechowski/meliphagid/analysis/2024-08-12/12-bwa/raw'
GENOME='/n/holylfs04/LABS/edwards_lab/Lab/smorzechowski/meliphagid/ReferenceAssemblies/Nleucotis_hifi_v1.0_hc_sm_fx_scaffolded_PAR_masked.fasta'
INDEXBASE='Nleucotis_hifi_v1.0_hc_sm_fx_scaffolded_PAR_masked'
FASTQ1=$1
FASTQ2=$2
CPU=16
FLOWCELL='LH00541'
LANE=$3
FLOWCELL_BARCODE='22KV7YLT3'
SAMPLE=$4
PLATFORM='Illumina'
LIBPREP='P536'
```

Next I activate my mamba environment for samtools and define filepaths for where I installed PICARD, BWA, and saved the perl script for getting alignment statistics.
```
module purge
module load python
source activate samtools
#source activate BamUtil

module load jdk/20.0.1-fasrc01

PICARD='/n/home09/smorzechowski/bin/picard/build/libs/picard.jar'
# this comes from ArimaGenomics mapping_pipeline
STATS='/n/home09/smorzechowski/bin/get_stats.pl'
BWA_PATH='/n/home09/smorzechowski/bin/bwa'
#BWA_PATH='/n/home09/smorzechowski/bin/bwa-mem2-2.2.1_x64-linux'
```

Next I check for directories I defined above and create them if they don't already exist. I also create the `bwa index` files for the reference genome if they don't already exist. Index files help BWA navigate and align short reads to the reference genome more efficiently.

Finally, I call BWA and pipe the output to `samtools` to create a `.bam` alignment file as BWA runs. It is best practice to automatically convert to `.bam` format rather than creating intermediate `.sam` files, which are very large and unwieldy. This BWA command also adds read group information defined in the input parameters above. It is important to include read group metadata such as: 
- FLOWCELL (from fastq file)
- LANE (from fastq file or your sequencing records)
- FLOWCELL_BARCODE (from fastq file)
- SAMPLE (from your sequencing records)
- PLATFORM (e.g. Illumina, PacBio, Oxford Nanopore)
- LIBPREP (from Bauer Core sequencing records)

The 'read group' metadata information gets associated with the mapped reads and is useful in a variety of contexts. For example, if you have reads that were sequenced on different lanes or with different library preps, it's important to identify which reads came from which lane to avoid technical artifacts that might occur if you combine all reads together without properly attributing what lane or library prep they came from.  

```
#Check output directories exist & create them as needed
[ -d $RAW_DIR ] || mkdir -p $RAW_DIR


# create bwa index if necessary - only do this once
[ -f ${INDEXBASE}.sa ] || $BWA_PATH/bwa index -p $INDEXBASE $GENOME

# call bwa and convert to bam, directly sort
# https://www.biostars.org/p/319730/
# -M : "mark shorter split hits as secondary
# -R : Read group information
# '@RG\tID:$ID\tSM:$SM\tLB:$LB\tPU:$PU\tPL:$PL'

$BWA_PATH/bwa mem $INDEXBASE $FASTQ1 $FASTQ2 -t $CPU -M -R $(echo "@RG\tID:$FLOWCELL"_"$LANE\tSM:$SAMPLE\tLB:$LIBPREP\tPU:$FLOWCELL_BARCODE"_"$LANE"_"$SAMPLE\tPL:$PLATFORM") \
 | samtools sort -@ $CPU --output-fmt BAM -o $RAW_DIR/${SAMPLE}_${LANE}_bwa_sort.bam
```

Finally I validate the bam files with PICARD to make sure that they are complete, not corrupted, formatted correctly:

```
# validate the bam files
java -Xmx30G \
-jar $PICARD ValidateSamFile \
      I=$RAW_DIR/${SAMPLE}_${LANE}_bwa_sort.bam \
      MODE=SUMMARY
```

## ANGSD pipeline

### Flagging duplicates and clipping reads


