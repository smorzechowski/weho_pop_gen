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


# $1 = R1 reads
# $2 = R2 reads

TRIMGALORE="/n/home09/smorzechowski/bin/TrimGalore-0.6.6/trim_galore"

# SEE OPTIONS $TRIMGALORE --help

# remove retained_unpaired, length=36 for paired
# set quality to 0 to include all base qualities, which ANGSD can handle and filter
# removing bases with quality less than 20 is standard for GATK and other variant callers

$TRIMGALORE --quality 20 --phred33 --fastqc --paired --output_dir trimmed_reads_all_qual --length 36 --stringency 3 -e 0.1 $1 $2
```

To run trimgalore on many individuals, I created a trimgalore.runscript using sed/awk to submit each trimgalore job to SLURM
See trimgalore.runscript in files.

```
```


