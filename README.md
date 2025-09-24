# Scripts associated with a population genomic analysis of White-eared Honeyeaters (*Nesoptilotis leucotis*)

This readme provides an overview of the genomic analyses conducted in this project focusing on genotype-environment associations in White-eared Honeyeaters to identify genomic regions involved in adaptation to aridity. I primarily ran the ANGSD pipeline for variant calling because I sequenced 72 individuals at an average depth of 6x coverage, which is in the low-to-mid coverage range suited for the genotype likelihood framework of ANGSD. Howevever, I include my scripts for hard-calling variants with GATK as well. I deployed a manual variant calling pipeline for the purpose of robustly analyzing variants on neo-sex chromosomes in males and females. 

I am quite excited by new automated Snakemake pipelines like [snpArcher](https://github.com/harvardinformatics/snpArcher) (implementing GATK) and [PopGLen](https://github.com/zjnolen/PopGLen) (implementing ANGSD). I believe snpArcher has an option for specifying ploidy, which is excellent, but I haven't yet implemented it or checked out the capabilities of PopGLen. I think routinely including sex chromosomes in population genomic and conservation genomic analyses would be a big step forward, especially if they could be fully incorporated into automated pipelines, which researchers with less familiarity with sex chromosomes could easily utilize and understand.   


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
- [local_pcangsd](https://github.com/alxsimon/local_pcangsd)
  - [lostruct](https://github.com/jguhlin/lostruct-py)
- [ngsRelate](https://github.com/ANGSD/NgsRelate)
- [ngsLD](https://github.com/fgvieira/ngsLD)
- [vcftools](https://vcftools.github.io/index.html)
- [bcftools](https://samtools.github.io/bcftools/bcftools.html)
- [qualimap](http://qualimap.conesalab.org/)
- [R (v4.3.2)](https://www.r-project.org/)


I was greatly assisted by helpful discussion with Dr. Elsie Shogren and her excellent [github repository](https://github.com/ehshogren/MyzomelaPopulationGenomics), the comprehensive [lsWGS tutorials](https://github.com/nt246/lcwgs-guide-tutorial) from the Therkildsen Lab, as well as a [custom python script](https://github.com/drewschield/Z-chromosome_analysis_hirundo/blob/main/scripts/identify_female_Zhet_sites.py) by Dr. Drew Schield, recommended to me  by Elsie. This script parses a VCF file from GATK to identify and remove spurious heterozygote variant calls on the Z chromosome in females. Elsie adapted this script to identify and remove spurious heterozygote calls on the W chromosome in females as well.       

I was also greatly facilitated by suggestions and helpful discussion with Dr. Teresa Pegan, who recommended a [github repository](https://github.com/alxsimon/local_pcangsd) that combines PCAngsd with lostruct in python to run local PCA in sliding windows across the genome to detect outliers in PCA space (often indicative of structural variants like inversions).


# Contents
- [Jobscripts and runscripts](#jobscripts-and-runscripts)
- [Order of operations for each pipeline](#order-of-operations-for-each-pipeline)
- [Genome assembly and curation](#genome-assembly-and-curation)
- [Adapter trimming and read mapping](#adapter-trimming-and-read-mapping)
- [ANGSD pipeline](#angsd-pipeline)
  - [Merging bam files](#merging-bam-files)
  - [Removing duplicates and clipping reads](#Removing-duplicates-and-clipping-reads)
  - [Indel realignment](#indel-realignment)
  - [Estimating coverage to verify sex chromosome complement](#estimating-coverage-to-verify-sex-chromosome-complement)
  - [Calculating genotype likelihoods](#calculating-genotype-likelihoods)
  - [PCA and population structure](#pca-and-population-structure)
  - [Local PCA with `local_pcangsd`](#local-pca-with-local_pcangsd)
  - [Creating imputed vcf files for LEA](#creating-imputed-vcf-files-for-lea)
  - [Genome-wide summary statistics](#genome-wide-summary-statistics)
- [GATK pipeline](#gatk-pipeline)
  - [Merging bam files and marking duplicates](#merging-bam-files-and-marking-duplicates)
  - [Split genome into intervals](#split-genome-into-intervals)
  - [Running haplotype caller](#running-haplotype-caller)
  - [Setting up genomicsDB for import](#setting-up-genomicsdb-for-import)
  - [Running genotypeGVCF](#running-genotypeGVCF)
- [Genotype-environment association analysis](#genotype-environment-association-analysis)
  - [sNMF analysis to estimate K ancestry proportions](#snmf-analysis-to-estimate-k-ancestry-proportions)
  - [LEA analysis](#lea-analysis)
- [Enrichment of candidate climate genes](#enrichment-of-candidate-climate-genes)
- [Association between body size and environmental variables](#association-between-body-size-and-environmental-variables)

## Jobscripts and runscripts

For most analyses, I provide both jobscripts and runscripts. The jobscript houses the self-contained code to run the analysis on the cluster, whereas the runscript is a simple bash script that provides the filepaths to the input files and submits the jobscript to the SLURM scheduler on the cluster with the `sbatch` command. I partition my analyses into these two scripts, thanks to a recommendation from Daren Card, because it facilitates record-keeping. The jobscript remains unchanged because I have used positional parameters (e.g. $1 $2) to make the jobscript universal. Parameter $1 refers to the first input provided to `sbatch`; parameter $2 refers to the second input, and so on. Each time I run the analysis, I add a new `sbatch` command to the runscript:

```
#!/bin/bash

input1='file/path/filename1.txt'
input2='file/path/filename2.txt'
input3='file/path/filename3.txt'

# Run X analysis again
sbatch X.jobscript $input1 $input2 $input3
```

## Order of operations for each pipeline

### GATK 
1. trimgalore.jobscript
2. bwa_sort.jobscript
3. merge_dedup_flag.jobscript
4. split_fastaindex.py
5. haplotypecaller.jobscript
6. genomicsDB_import.jobscript
7. genotypeGVCF.jobscript
8. Combine VCF files
9. Filter VCF files
10. Run PCA, etc.

### ANGSD 
1. trimgalore.jobscript
2. bwa_sort.jobscript
3. merge_dedup_clip.jobscript
4. indel_realign.jobscript


## Genome assembly and curation

I generated HiFi PacBio reads on a Revio SMRT Cell from a female White-eared Honeyeater to create a long-read reference genome containing both the neo-Z and neo-W sex chromosomes. I used Heng Li's fantastic assembler `hifiasm` to conduct *de novo* diploid assembly of the HiFi reads. I used a handy snakemake pipeline from Danielle Khost of the Harvard Informatics group (https://github.com/harvardinformatics/pacbio_hifi_assembly) to run `hifiasm` on the cluster. 

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

Hifiasm is able to produce two haplotype assemblies (hap1 and hap2) as well as a primary haploid assembly that collapses the bubbles in the graph. Since I am not doing pangenome analysis, I can use the primary assembly to do short-read alignment for variant calling. However, it is often necessary to manually curate the sex chromosomes, especially in a species like White-eared Honeyeater, which contains neo-sex chromosomes. I found that the primary assembly did not contain the full phased sequence for the neo-Z or neo-W because the algorithm incorrectly purged haplotigs that it characterized as redundant but were actually minimally diverged regions of the neo-sex chromosomes in this species. It is not surprising that it purged the redundant contig(s) containing the new PAR because it is not diverged at all between the neo-sex chromosomes. However, it also purged regions where recombination suppression has occurred on the neo-W, leading to nascent divergence between the Z and W gametologs. 

To make sure that all W-linked and Z-linked contigs were included in the final assembly, I had to identify the Z- and W-linked contigs in all assemblies with a combination of sex differences in coverage and synteny analysis.

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
Here a provide a comprehensive summary of all the scripts I generated to run ANGSD and estimate genome summary statistics within genotype likelihood framework.

### Merging bam files

All my 72 lcWGS samples were sequenced across two lanes of an Illumina NovaSeq X in order to have the coverage I requested (~6x for autosomes per individual). I decided not to go much lower than that because I wanted decent coverage (>3x) on the Z and W chromosomes in females. Also, sometimes individuals have to be excluded if they are too low relative to others due to the normal stochasticity that happens during library prep and sequencing. So, after I mapped and validated the bam files, I had to merge the two bam files per individual:


```
# check output directories & create them as needed
[ -d $RAW_DIR ] || mkdir -p $RAW_DIR
[ -d $MERGE_DIR ] || mkdir -p $MERGE_DIR

# merge the two bam files for each sample
samtools merge -r $MERGE_DIR/${SAMPLE}_bwa_merge.bam $RAW_DIR/${SAMPLE}_L003_bwa_sort.bam $RAW_DIR/${SAMPLE}_L004_bwa_sort.bam

# view the header to look at read groups
#samtools view -H $MERGE_DIR/${SAMPLE}_bwa_merge.bam

# sort the merged file again
samtools sort $MERGE_DIR/${SAMPLE}_bwa_merge.bam -@ $CPU --output-fmt BAM -o $MERGE_DIR/${SAMPLE}_bwa_merge_sort.bam

# index the sorted merged bam files
# not sure if this is necessary?
#samtools index $RAW_DIR/${SAMPLE}_bwa_merge_sort.bam

# get summary of the alignment
samtools flagstat $MERGE_DIR/${SAMPLE}_bwa_merge_sort.bam > $MERGE_DIR/${SAMPLE}_bwa_merge_sort.bam.flagstat

```


### Removing duplicates and clipping reads

Next, it is necessary to flag/remove PCR and optical duplicates, which are a normal part of the Illumina short-read sequencing process. I used Picard to mark and *remove* the duplicates (`REMOVE_DUPLICATES=true`) because I was not certain if ANGSD interprets flags in the same way that GATK does. For GATK, it is sufficient to mark reads as being duplicates, rather than removing them; once flagged, GATK will not include them in the variant calling process. Upon further digging, it seems that either flagging or removing duplicates works in ANGSD as well as GATK, I just went with the conservative option of removing them altogether.

```
# Mark duplicates with Picard
#-XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/ \
java -Xmx30G \
-jar $PICARD MarkDuplicates \
TMP_DIR=$TMP_DIR \
INPUT=$MERGE_DIR/${SAMPLE}_bwa_merge_sort.bam \
OUTPUT=$DEDUP_DIR/${SAMPLE}_bwa_merge_dedup_flag.bam \
METRICS_FILE=$DEDUP_DIR/${SAMPLE}_bwa_merge_dedup_metrics.txt \
REMOVE_DUPLICATES=true \
TAGGING_POLICY=All \
ASSUME_SORTED=true \

```

Next, I needed to clip overlapping reads for downstream ANGSD analyses because ANGSD is not able to account for the fact that forward and reverse reads may be overlapping unless they are explicitly flagged as such. ANGSD tutorials generally recommend `bamUtil clipOverlap` for this, but I was finding that running this software created a lot of invalid CIGAR strings, etc. when I ran `ValidateSamFile`. It may have been because of software conflict, I'm not sure. [Others](https://github.com/statgen/bamUtil/issues/72) have discovered this issue as well. Regardless, I found another way to clip overlapping reads with `fgbio` [ClipBam](https://fulcrumgenomics.github.io/fgbio/tools/latest/ClipBam.html), which seems to have worked great! 

Note: I chose to softclip the overlaps rather than hardclipping. From what I understand, ANGSD can read the CIGAR flags and ignore overlapping softclipped regions accordingly. A more conservative option is hardclipping all overlaps, however, I was trying to replicate the behavior of `bamUtil clipOverlap`, which to the best of my knowledge and digging, does do softclipping (see comment from mktrost [here](https://github.com/statgen/bamUtil/issues/15) for instance). However, if anyone finds any information to the contrary, I am happy to stand corrected!

```

# need to sort reads by query name for fgbio ClipBam - see user manual
# remove -u which uncompresses, I don't think we want that
samtools sort -n -@ $CPU $DEDUP_DIR/${SAMPLE}_bwa_dedup.bam -o $DEDUP_DIR/${SAMPLE}_bwa_dedup_namesort.bam

# softclip reads using ClipBam from fgbio
fgbio -Xmx45G ClipBam \
    -i $DEDUP_DIR/${SAMPLE}_bwa_dedup_namesort.bam \
    -o $CLIP_DIR/${SAMPLE}_bwa_dedup_clip.bam \
    -m $CLIP_DIR/${SAMPLE}_clipbam.metrics \
    -r $GENOME \
    --clip-overlapping-reads=true \
    -c Soft

# Sort the reads normally again
samtools sort -@ $CPU $CLIP_DIR/${SAMPLE}_bwa_dedup_clip.bam -o $CLIP_DIR/${SAMPLE}_bwa_dedup_clip_sort.bam

```

Finally, I ran Picard's `FixMate` to make sure that mate-pair information such as `mate position`, `insert size`, `mate strand orientation` and `proper pair flag` was consistent across paired reads. Then, I reindexed, ran output stats, and validated the bam files one last time.

```
PICARD='/n/home09/smorzechowski/bin/picard/build/libs/picard.jar'
# this comes from ArimaGenomics mapping_pipeline
STATS='/n/home09/smorzechowski/bin/get_stats.pl'

[ -d $CLIP_DIR ] || mkdir -p $CLIP_DIR

# fix mate information
# output will be sorted automatically the same as input, e.g. by coordinate
java -Xmx30G \
-jar $PICARD FixMateInformation \
       I=$CLIP_DIR/${SAMPLE}_bwa_dedup_clip_sort.bam \
       O=$CLIP_DIR/${SAMPLE}_bwa_dedup_clip_sort_fixmate.bam \
       ADD_MATE_CIGAR=true

# Index the final bam files again
samtools index $CLIP_DIR/${SAMPLE}_bwa_dedup_clip_sort_fixmate.bam

# output stats on mapping using ArimaGenomics perl script
perl $STATS $CLIP_DIR/${SAMPLE}_bwa_dedup_clip_sort_fixmate.bam > $CLIP_DIR/${SAMPLE}_bwa_dedup_clip_sort_fixmate.bam.stats


# validate the fixed mate bam files
java -Xmx30G \
-jar $PICARD ValidateSamFile \
      I=$CLIP_DIR/${SAMPLE}_bwa_dedup_clip_sort_fixmate.bam \
      MODE=SUMMARY
```

### Indel realignment

Next, following the excellent tutorial on analyzing lcWGS data from the Therkildsen Lab (see https://github.com/nt246/lcwgs-guide-tutorial), I performed indel realignment. Often short reads are mismapped in the vicinity of indels, which can lead to spurious variant calling. Thus, the GATK tools `RealignerTargetCreater` which identifies all indel regions, and `IndelRealigner` which does the realignment, are dedicated to fixing mismapping around indels. This is important for ANGSD, however, these tools are no longer supported in GATK > 3.7, thus you have to go to the archives to download and use the correct version: `GATK 3.7`. 

```
## Create dict file for reference genome
#module load jdk/20.0.1-fasrc01
#PICARD='/n/home09/smorzechowski/bin/picard/build/libs/picard.jar'
#java -jar $PICARD CreateSequenceDictionary \
#      R=$GENOME \
#      O=/n/holylfs04/LABS/edwards_lab/Lab/smorzechowski/meliphagid/ReferenceAssemblies/Nleucotis_hifi_v1.0_hc_sm_fx_scaffolded_PAR_masked.dict

## Realign around indels
# This is done across all samples at once

## Create list of potential indels
# this can be multithreaded, but IndelRealigner cannot

java -Xmx60G -jar /n/home09/smorzechowski/.conda/envs/gatk_3.7/opt/gatk-3.7/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R $GENOME \
-I $BAMLIST \
-o $OUTDIR'/all_samples_for_indel_realigner_incl_2final.intervals' \
-drf BadMate \
-nt 20


## Run the indel realigner tool
java -Xmx25G -jar /n/home09/smorzechowski/.conda/envs/gatk_3.7/opt/gatk-3.7/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R $GENOME \
-I $BAMDIR/${sample}_bwa_dedup_clip_sort_fixmate.bam \
-targetIntervals $OUTDIR'/all_samples_for_indel_realigner_incl_2final.intervals' \
--consensusDeterminationModel USE_READS  \
--nWayOut _realigned_final.bam
```

### Estimating coverage to verify sex chromosome complement

Next, before you actually calculate genotype likelihoods, it is important to verify the sex chromosome complement of all individuals you sequenced. Misclassification can happen, especially with species that are monomorphic or monochromatic like White-eared Honeyeaters. Separating males (homogametic ZZ) and females (heterogametic ZW) is particularly important when doing variant calling on the sex chromosomes because males are diploid for Z and females are haploid for both Z and W. This will lead to spurious variant calling if you combine both sexes and try to call variants on the Z, for instance. What I was suprised to find out, however, is that ANGSD does not appear to be able to accomodate variant calling on haploid Z or W chromosomes, which is quite unfortunate. This meant that I had to restrict my samples to males when conducting analyses involving the Z chromosome.

An easy way to verify the sex chromosome complement is to estimate coverage of reads mapping to autosomes, Z, and W chromosomes for each individual. If autosomal coverage = Z coverage, I characterize this as a ZZ male individual. If autosomal coverage = 2 x Z coverage, this is characterized as a ZW individual because Z coverage is half that of autosomes. Similarly if W coverage is around zero, that is also characterized as a ZZ male individual. Finally, if autosomal coverage = 2 x W coverage, that is also characterized as a ZW female individual, the same as you would expect with Z coverage in females. 

It is worth noting that sex chromosome complement does not necessarily always align with other sex variables, depending on the individual and the species. However, for the purpose of this work on sex chromosomes, I usually talk about sex in the context of the sex chromosome complement. 

I used qualimap to estimate read depth.

```
module load python
source activate /n/holylabs/LABS/edwards_lab/Users/smorzechowski/conda/qualimap

sample=$1

bam_path=/n/netscratch/edwards_lab/Lab/smorzechowski/meliphagid/analysis/2024-11-03/06-bwa/clip
out_path=/n/netscratch/edwards_lab/Lab/smorzechowski/meliphagid/analysis/2024-11-03/10-qualimap/qualimap_final

qualimap bamqc -bam ${bam_path}/${sample}_bwa_dedup_clip_sort_fixmate.bam -sd -c -nw 400 -hm 3 -outdir $out_path/${sample} --java-mem-size=25G

```

Then, I extracted the coverage information using an awk command. 

```
#!/bin/bash/

sample_names='samplenames_all74.txt'

cat $sample_names |
while read sample; do
awk '/>>>>>>> Coverage per contig/{flag=1;next}/ /{flag=0}flag' ./qualimap_final/${sample}/genome_results.txt | \
awk '{print $1,$2,$3,$4,$5}' | awk '$(NF+1)= "'$sample'"' > ${sample}_coverage_results_final.txt
done

cd results_final
cat *final.txt > all_coverage_results_final.txt
awk 'NF >1' all_coverage_results_final.txt > all_coverage_results_cleaned_final.txt

```

Finally, in R, I plotted the coverage per individual for an autosome (Chr 1), and the W chromosome.

```
# all_coverage_results_cleaned_final.txt contains the two individuals from Moonbi sequenced seperately

data <- read.table("all_coverage_results_cleaned_final.txt",header=F,sep=" ",comment.char = "")

data_filt <- data[data$V1=="scaffold_5_RagTag" |
                  data$V1=="scaffold_2_RagTag" |
                  data$V1=="scaffold_1_RagTag",]


autos <- data[!is.element(data$V1,c("scaffold_2_RagTag","scaffold_5_RagTag")),]

auto_sum <- autos %>% group_by(V6)%>%
  summarise(mean_coverage = mean(V4))

mean(auto_sum$mean_coverage)
max(auto_sum$mean_coverage)
min(auto_sum$mean_coverage)


auto1a <- data[data$V1=="scaffold_1_RagTag",]
auto1a_sum <- auto1a %>% group_by(V6)%>%
  summarise(mean_coverage = mean(V4))

mean(auto1a_sum$mean_coverage)
max(auto1a_sum$mean_coverage)
min(auto1a_sum$mean_coverage)



ggplot(data_filt,aes(V6,V4))+
  geom_bar(stat='identity')+
  facet_grid(.~V1)

# plot the W chromosome coverage
data_W <- data[data$V1=="scaffold_5_RagTag",]
ggplot(data_W,aes(V6,V4))+
  geom_bar(stat='identity')+
  coord_flip()+
  ylab('Neo-W Coverage (excluding new PAR)')+
  xlab("")

# plot the autosome coverage
data_1a <- data[data$V1=="scaffold_1_RagTag",]
ggplot(data_1a,aes(V6,V4))+
  geom_bar(stat='identity')+
  coord_flip()+
  ylab('1A coverage (autosomal)')+
  xlab("")

```


### Calculating genotype likelihoods

I calculated genotype likelihoods using the program ANGSD. I did this on the full dataset (all individuals, autosomes only) and the male dataset (autosomes + neo-Z chromosome). I employed various filters and quality cut-offs and produced genotype likelihoods in the beagle file format. 

Example: Calculating genotype likelihoods for autosomes and the full dataset:

```
BASEDIR='/n/netscratch/edwards_lab/Lab/smorzechowski/meliphagid/analysis/2024-11-03/14-angsd'
ANGSD='/n/home09/smorzechowski/bin/angsd/angsd'
REF='/n/holylfs04/LABS/edwards_lab/Lab/smorzechowski/meliphagid/ReferenceAssemblies/Nleucotis/Nleucotis_hifi_v1.0_hc_sm_fx_scaffolded_PAR_masked.fasta'
RM_DIR='/n/holylfs04/LABS/edwards_lab/Lab/smorzechowski/meliphagid/ReferenceAssemblies/Nleucotis/Nleu_final_genome_softmask'

PCANGSD='/n/home09/smorzechowski/bin/pcangsd/pcangsd/pcangsd.py'

# index the sites file, do this just once
$ANGSD sites index $RM_DIR'/Nleucotis_hifi_v1.0_hc_sm_fx_scaffolded_PAR_masked_sites_keep_1based.txt'

bamlist=$1
region=$2
out=$3

# activate dependencies for pcangsd
module load python
source activate pcangsd
module purge

# PCA based on all SNPs without LD pruning first
# minor allele freq > 0.05

$ANGSD -b $BASEDIR'/sample_lists/'$bamlist \
-anc $REF \
-rf $BASEDIR'/regions/'$region \
-out $BASEDIR'/results/'$out \
-minMapQ 30 \
-minQ 30 \
-doMaf 1 \
-minMaf 0.05 \
-SNP_pval 2e-6 \
-GL 1 \
-doGlf 2 \
-doMajorMinor 1 \
-skipTriallelic 1 \
-doPost 1 \
-doIBS 1 \
-doCounts 1 \
-doCov 1 \
-makeMatrix 1 \
-P 8 \
-sites $RM_DIR'/Nleucotis_hifi_v1.0_hc_sm_fx_scaffolded_PAR_masked_sites_keep_1based.txt' \
-remove_bads 1 \
-uniqueOnly 1 \

```

Then, I ran PCAngsd to calculate PCA and visualize population structure. 

```
# Run PCAngsd
pcangsd \
-b $BASEDIR'/results/'$out'.beagle.gz' \
-o $BASEDIR'/results/'$out'_pcangsd_thinned' \
--geno 0.1 \
--threads 8 \

```

### PCA and population structure

### Local PCA with `local_pcangsd`

The python script I adapted from [Alexis Simon](https://www.normalesup.org/~asimon/projects/local_pcangsd.html) to combine PCAngsd and lostruct can be found [here](https://github.com/smorzechowski/weho_pop_gen/blob/master/local_pcangsd_scripts/local_pcangsd_script.py).

Basically, the python script runs PCAngsd separately in pre-defined windows (e.g. a window size of 50,000 bp with a minimum number of 500 variants, for example) across each scaffold/chromosome. It calculates the pair-wise distance between windows, summarizes the variation with Principal Coordinates Analysis (PCoA, which is equivalent to classical metric multidimensional scaling, or MDS), and plots the first two axes (MDS1 and MDS2) in bivariate space. It finds outlier regions where the MDS scores are different from the rest of the chromosome/scaffold. It then merges those outlier windows, runs PCAngsd on this region, and plots the first two PCA coordinates to visualize the population structure in this outlier region. A PCA with three distinct clusters defined across PC1 and no differentiation across PC2 is a potential indication of an inversion polymorphism.  

This program requires a beagle file as input. I split the beagle files by chromosomes/scaffold and ran local_pcangsd separately on each. The input is compressed in zarr format, and the results are written to a zarr file. 

```
cat unplaced_scaffolds_beagles.txt | while read line
do
prefix=`echo $line |sed 's/.beagle.gz//g'`
input="/n/netscratch/edwards_lab/Lab/smorzechowski/meliphagid/analysis/2024-11-03/14-angsd/results/$line"
store="/n/netscratch/edwards_lab/Lab/smorzechowski/meliphagid/analysis/2024-11-03/16-local_pcangsd/zarr_dir/${prefix}.zarr"
zarr="/n/netscratch/edwards_lab/Lab/smorzechowski/meliphagid/analysis/2024-11-03/16-local_pcangsd/results/${prefix}_results.zarr"
tmp="/n/home09/smorzechowski/tmp_local_pcangsd"
outdir="/n/netscratch/edwards_lab/Lab/smorzechowski/meliphagid/analysis/2024-11-03/16-local_pcangsd/plots/"
col=0
pop='/n/netscratch/edwards_lab/Lab/smorzechowski/meliphagid/analysis/2024-11-03/16-local_pcangsd/all_samples_pops.tsv'
wind=10000
var=50

python /n/home09/smorzechowski/local_pcangsd_script.py $input $store $zarr $tmp $prefix $outdir $col $pop $wind $var
done
```

### Creating imputed vcf files for LEA

First, I ran ANGSD with all the necessary filters (excluding sites with more than 15% missing data, for example, using `-minInd 55` for the full dataset and `-minInd 36` for the male dataset) the `-doBcf 1` option to create the .bcf file format needed for converting to .vcf.gz for the imputation program BEAGLE.

```
$ANGSD -b $BASEDIR'/sample_lists/'$bamlist \
-ref $REF \
-rf $BASEDIR'/regions/'$region \
-out $BASEDIR'/results/'$out \
-minMapQ 30 \
-minQ 30 \
-doMaf 1 \
-minMaf 0.05 \
-SNP_pval 1e-6 \
-GL 1 \
-doGlf 2 \
-doMajorMinor 1 \
-skipTriallelic 1 \
-doPost 1 \
-doIBS 1 \
-doCounts 1 \
-doCov 1 \
-makeMatrix 1 \
-P 8 \
-sites $RM_DIR'/Nleucotis_hifi_v1.0_hc_sm_fx_scaffolded_PAR_masked_sites_keep_1based.txt' \
-remove_bads 1 \
-uniqueOnly 1 \
-minInd 36 \
-doBcf 1 \
```

I filtered for biallelic SNPs and converted the bcf format to vcf.gz format with bcftools and then filtered for the desired maximum mean depth. 

```
module load python
source activate bcftools
module purge

#input=$1
#input='/n/netscratch/edwards_lab/Lab/smorzechowski/meliphagid/analysis/2024-11-03/14-angsd/results/Nleu_autos_lea.bcf'
input='/n/netscratch/edwards_lab/Lab/smorzechowski/meliphagid/analysis/2024-11-03/14-angsd/results/Nleu_autos_neoZ_lea_males.bcf'
prefix=$(echo "$input" | sed 's:.*/::; s/.bcf//g')

#prefix=$(basename "$input" .bcf)
#prefix=`echo $input |sed 's/.bcf//g'`

# only include snps and biallelic sites
bcftools view -v snps --max-alleles 2 --threads $SLURM_CPUS_PER_TASK $input -Oz -o $prefix'.vcf.gz'

module load python
source activate vcftools
module purge


#bgzip $prefix'.vcf'
tabix -p vcf $prefix'.vcf.gz'

vcf=$prefix'.vcf.gz'

vcftools --gzvcf $vcf --max-meanDP 9 --recode --stdout | bgzip > $prefix'_depth_filt.vcf.gz'

```

Then, with BEAGLE, I imputed variants for autosomes (all individuals) and autosomes plus the neo-Z chromosome (males only).

```
bin='/n/home09/smorzechowski/bin'
module load jdk/20.0.1-fasrc01

input=$1
#map=$2
out=$2

java -Xss5m -Xmx80g -jar $bin/beagle.27Jan18.7e1.jar \
gl=$input \
out=$out \
nthreads=16 \
```


Run script: 
```
#!/bin/bash

# Run for autosomes
#input='Nleu_autos_lea_depth_filt.vcf.gz'
#map='Nleu_autos_lea_depth_filt_plink_v2.map'
#output='Nleu_autos_lea_depth_filt_imputed'
#sbatch beagle_imputation.jobscript $input $output

# Run for males: autos and neoZ
input='Nleu_autos_neoZ_lea_males_depth_filt.vcf.gz'
output='Nleu_autos_neoZ_lea_males_depth_filt_imputed'
sbatch beagle_imputation.jobscript $input $output

```
### Genome-wide summary statistics

## GATK pipeline

I provide a summary of the GATK pipeline below, even though I did not end up using the vcf files with hard-called variants from GATK because the genotype likelihood framework from ANGSD was better suited to my low-to-mid coverage data.

### Marking duplicates

## Genotype-environment association analysis

### sNMF analysis to estimate K ancestry proportions

I used the `snmf()` function in the R package `LEA` to determine the best value of K ancestry proportions to use as latent factors in subsequent GEA analyses. 

```
# Calculate the cross entropy criterion to pick the best value of K ancestral populations

setwd('/n/netscratch/edwards_lab/Lab/smorzechowski/meliphagid/analysis/2024-11-03/19-lea')


library(LEA)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(stringr)
library(cowplot)
library(gridExtra)



gen.imp <- read.geno("Nleu_autos_lea_depth_filt_imputed_thin_plink.geno")

dim(gen.imp)

################
# running snmf #
################

project.snmf = snmf("Nleu_autos_lea_depth_filt_imputed_thin_plink.geno",
                    K = 1:10, 
                    entropy = TRUE, 
                    repetitions = 10,
                    project = "new")

# using project = "new" saves it to a file with .snmfProject extension

png("sNMF_K1_to_10.png", width = 6, height = 4, units = "in", res = 350)
# plot cross-entropy criterion of all runs of the project
plot(project.snmf, cex = 1.2, col = "lightblue", pch = 19)
dev.off()

```


### LEA analysis

For running latent factor mixed models in the R package `LEA` using the `lfmm2()` function, I first had to convert the input files to .geno format (`ped2geno`). I also converted the text file of environmental predictors to .env format with `write.env()`.

The environmental data can be summarized into the first two principal components, otherwise specific, relevant variables can be selected from the BioClim database to test for genotype-environment associations. 

The most basic LEA analysis consists of reading in the data, running the models with specified values of K (the number of latent factors), computing the multivariate tests of significance and adjusting the p-values with the BH procedure. 


```

gen.imp <- read.geno("Nleu_autos_lea_depth_filt_imputed_thin_plink.geno")

dim(gen.imp)

#pred <- read.env("Nleu_bioclim_variables_65ind_autos_ordered_noheader.env")
pred <- read.env("Nleu_bioclim_PC1_PC2_65ind_autos.env")

# Run LEA models
mod2 <- lfmm2(input = gen.imp, env = pred, K = 2)
#mod3 <- lfmm2(input = gen.imp, env = pred, K = 3)
#mod4 <- lfmm2(input = gen.imp, env = pred, K = 4)

# Run multivariate tests, accounting for shared variance (full covariance matrix)
# K=2
pv2.full <- lfmm2.test(object = mod2,
                 input = gen.imp,
                 env = pred,
                 full=TRUE)

# Calculate adjusted pvalues for the full model
pv2.full.q.values <- p.adjust(pv2.full$pvalues, method = "BH")
#pv3.full.q.values <- p.adjust(pv3.full$pvalues, method = "BH")

```

The most basic plotting of the results can be done in base R. 

```

# Target regions across the genome: chromosomes with putative inversions

Chr4_target = seq(from = 474812, to = 583562,by=1)
Chr6_target = seq(from = 583563, to = 650966,by=1)
Chr8_target = seq(from = 700446, to = 757547,by=1)
Chr9_target = seq(from = 757547, to = 773605,by=1)
Chr10_target = seq(from = 139148, to = 158627,by=1)
Chr13_target = seq(from = 172691, to = 201522,by=1)
Chr14_target = seq(from = 201523, to = 220168,by=1)
Chr15_target = seq(from = 220169, to = 230157,by=1)
Chr16_target = seq(from = 230158, to = 247440,by=1)
Chr17_target = seq(from = 247441, to = 263872,by=1)
Chr18_target = seq(from = 263872, to = 272643,by=1)
Con15l_target = seq(from = 773606, to = 789999,by=1)
Con45l_target = seq(from = 790000, to = 807105,by=1)

# Plot of all variables combined with BH correction K=2
plot(-log10(pv2.full.q.values), col = "grey", cex = .4, pch = 19)
abline(h = -log10(0.0000005), lty = 2, col = "darkred")

points(Chr4_target, -log10(pv2.full.q.values[Chr4_target]), col = "red")
points(Chr6_target, -log10(pv2.full.q.values[Chr6_target]), col = "blue")
points(Chr8_target, -log10(pv2.full.q.values[Chr8_target]), col = "lightblue")
points(Chr9_target, -log10(pv2.full.q.values[Chr9_target]), col = "magenta")
points(Chr10_target, -log10(pv2.full.q.values[Chr10_target]), col = "orange")
points(Chr13_target, -log10(pv2.full.q.values[Chr13_target]), col = "green")
points(Chr14_target, -log10(pv2.full.q.values[Chr14_target]), col = "yellow")
points(Chr15_target, -log10(pv2.full.q.values[Chr15_target]), col = "black")
points(Chr16_target, -log10(pv2.full.q.values[Chr16_target]), col = "darkgreen")
points(Chr17_target, -log10(pv2.full.q.values[Chr17_target]), col = "pink")
points(Chr18_target, -log10(pv2.full.q.values[Chr18_target]), col = "brown")
points(Con15l_target, -log10(pv2.full.q.values[Con15l_target]), col = "darkblue")
points(Con45l_target, -log10(pv2.full.q.values[Con45l_target]), col = "darkorange")


# Add a legend to the plot
legend("topleft", 
       legend = c("Chr4", "Chr6", "Chr8", "Chr9", "Chr10", 
                  "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", 
                  "Chr18", "Con15l", "Con45l"),
       col = c("red", "blue", "lightblue", "magenta", "orange", 
               "green", "yellow", "black", "darkgreen", "pink", 
               "brown", "darkblue", "darkorange"),
       pch = 16, # Use circles as point markers
       bty = "n") # Remove box around legend

```

## Enrichment of candidate climate genes

I collated a list of avian candidate climate genes from an extensive literature search. I also used a list of vertebrate candidate climate genes assembled in [Wollenberg Valero et al. 2022](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2656.13617). 

I used the R package `regioneR` to conduct enrichment tests of these candidate climate genes in putative inversions. The goal was to determine if these candidates are significantly more likely to be found in inversions compared to random expectations.

```
library(regioneR)
library(GenomicRanges)  # For creating GRanges objects
library(rtracklayer)
library(ggplot2)
library(patchwork)
library(cowplot)

setwd("~/PhD research/Neo sex chromosome/WEHE pop gen chapter/WEHE pop gen/regioneR")

gr_PCR <- import("inversions_updated.bed",format="BED")

gr_avian <- import("Nleu_avian_candidate_climate_genes_sorted_IDs.bed")
gr_vertebrate <- import("Nleu_vertebrate_pEAGs_genes_sorted_IDs.bed")

genome <- read.table("Nleucotis_hifi_v1.0_hc_sm_fx_scaffolded_PAR_masked.fasta.genome",header=F,comment.char = "")
genome$V2 <- as.numeric(genome$V2)

genome_vector <- setNames(genome$V2,genome$V1)
str(genome_vector)
genome_gr <- GRanges(
  seqnames = names(genome_vector),
  ranges = IRanges(start = 1, end = as.numeric(genome_vector))
)



pt1 <- overlapPermTest(
  A = gr_avian,         # Candidate regions
  B = gr_PCR,          # Annotation regions (e.g., inversions)
  genome = genome_gr, # Specify the genome (or you can provide a custom genome definition)
  ntimes = 1000             # Number of permutations to build the null distribution
)

null_mean <- mean(pt1$numOverlaps$permuted)
fold_enrichment <- 98 / null_mean

fold_enrichment

data <- data.frame(pt1$numOverlaps$permuted)

pt1

plot1 <- ggplot(data,aes(pt1.numOverlaps.permuted))+
  geom_density(trim=TRUE)+
  xlab("Overlap of avian candidate climate genes in MORs")+
  ylab("Density")+
  ggtitle("A.")+
#  geom_vline(xintercept=79.337,linetype="dashed",color='black')+
  geom_vline(xintercept=98,color='red',linewidth=2)+
  theme_minimal()+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=14),
        title=element_text(size=20))
  

plot1



pt2 <- overlapPermTest(
  A = gr_vertebrate,         # Candidate regions
  B = gr_PCR,          # Annotation regions (e.g., inversions)
  genome = genome_gr, # Specify the genome (or you can provide a custom genome definition)
  ntimes = 1000             # Number of permutations to build the null distribution
)

head(pt2$numOverlaps)


null_mean <- mean(pt2$numOverlaps$permuted)
fold_enrichment <- 143 / null_mean

pt2
fold_enrichment

data <- data.frame(pt2$numOverlaps$permuted)

plot2 <- ggplot(data,aes(pt2.numOverlaps.permuted))+
  geom_density(trim=TRUE)+
  xlab("Overlap of vertebrate candidate climate genes in MORs")+
  ylab("Density")+
  ggtitle("B.")+
 # geom_vline(xintercept=95.497,linetype="dashed",color='black')+
  geom_vline(xintercept=143,color='red',linewidth=2)+
  theme_minimal()+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=14),
        title=element_text(size=20))

final_plot <- plot1/plot2
final_plot



```

## Association between body size and environmental variables


