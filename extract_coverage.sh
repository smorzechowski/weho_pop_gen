#!/bin/bash/

#sample_names='samplenames.txt'
#sample_names='samplenames_test.txt'
#sample_names='samplenames_remaining.txt'
sample_names='samplenames_all74.txt'

cat $sample_names |
while read sample; do
awk '/>>>>>>> Coverage per contig/{flag=1;next}/ /{flag=0}flag' ./qualimap_final/${sample}/genome_results.txt | \
awk '{print $1,$2,$3,$4,$5}' | awk '$(NF+1)= "'$sample'"' > ${sample}_coverage_results_final.txt
done

# cd results_final
# cat *final.txt > all_coverage_results_final.txt
# awk 'NF >1' all_coverage_results_final.txt > all_coverage_results_cleaned_final.txt

