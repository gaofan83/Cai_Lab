#!/bin/bash

export PATH=/home/fgao/software/STAR-2.7.8a/bin/Linux_x86_64/:$PATH

# star alignment
while read line;
do
ext1="_R1.fastq.gz"
ext2="_R2.fastq.gz"
file1=$line$ext1
file2=$line$ext2

STAR --runThreadN 32 \
--genomeDir /home/fgao/reference/STARIndex_GRCm38/ \
--readFilesIn ./$file1 ./$file2 \
--readFilesCommand zcat \
--outFileNamePrefix $line \
--limitBAMsortRAM 50000000000 \
–outFilterMultimapNmax 2 \
--outSAMstrandField intronMotif \
--outSAMtype BAM SortedByCoordinate

samtools index ${line}Aligned.sortedByCoord.out.bam
done < sample_ID1.txt
