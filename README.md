# Order of operations

## GATK pipeline
1. trimgalore.jobscript
2. bwa_sort.jobscript
3. merge_dedup_flag.jobscript
4. haplotypecaller.jobscript
5. genomicsDB_import.jobscript
6. genotypeGVCF.jobscript
7. Combine VCF files
8. Filter VCF files
9. Run PCA

## ANGSD pipeline
1. trimgalore.jobscript
2. bwa_sort.jobscript
3. merge_dedup_clip.jobscript
4. indel_realign.jobscript
