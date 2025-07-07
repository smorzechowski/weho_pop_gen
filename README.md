# Scripts associated with a population genomic analysis of White-eared Honeyeaters (*Nesoptilotis leucotis*)

This readme provides an overview of the genomic analyses conducted in this project. I primarily ran the ANGSD pipeline for variant calling because I sequenced 72 individuals at an average depth of 6x coverage, which is in the low-to-mid coverage range better suited for the genotype likelihood framework of ANGSD. Howevever, I include the manual GATK pipeline for hard calling variants as well.  

The software and programs used in this project include:

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
- [RepeatModeler](https://github.com/Dfam-consortium/RepeatModeler/tree/master)
- [RepeatMasker](https://github.com/Dfam-consortium/RepeatMasker)
- [R (v4.3.2)](https://www.r-project.org/)


# Contents
- [Order of opererations for each pipeline](#order-of-operations-for-each-pipeline)
- [Adapter trimming and read mapping](#adapter-trimming-and-read-mapping)


## Order of operations for each pipeline

### GATK pipeline
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

### ANGSD pipeline
1. trimgalore.jobscript
2. bwa_sort.jobscript
3. merge_dedup_clip.jobscript
4. indel_realign.jobscript


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

Finally I validate the bam files with PICARD to make sure that they are complete, not corrupted, formatted correctly. 
```
# validate the bam files
java -Xmx30G \
-jar $PICARD ValidateSamFile \
      I=$RAW_DIR/${SAMPLE}_${LANE}_bwa_sort.bam \
      MODE=SUMMARY

```



