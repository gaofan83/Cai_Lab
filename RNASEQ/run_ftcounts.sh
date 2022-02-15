#!/bin/bash

mkdir -p ./ft_out
mkdir -p ./ft_seq
while read line;
do
ext="Aligned.sortedByCoord.out.bam"
file_out=$line$ext
~/software/subread-1.6.4-Linux-x86_64/bin/featureCounts -T 48 -t exon -p -g gene_id -a /home/fgao/refdata-cellranger-mm10-3.0.0/genes/genes.gtf -o ./ft_out/ft_$line.txt $file_out
done < sample_ID1.txt

#-s 1 stranded; -s 2 reverse stranded;

while read line;
do
cat ./ft_out/ft_$line.txt | awk '{if($1!="#" && $1!="Geneid")print $1"\t"$7}' > ./ft_seq/ftcounts_$line.txt
done < sample_ID1.txt
