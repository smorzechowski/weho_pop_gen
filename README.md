# Order of operations

## GATK pipeline
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

## ANGSD pipeline
1. trimgalore.jobscript
2. bwa_sort.jobscript
3. merge_dedup_clip.jobscript
4. indel_realign.jobscript
