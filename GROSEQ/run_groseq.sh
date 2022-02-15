cutadapt -l 32 -m 32 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o rep1_trim.fastq.gz rep1.fastq.gz

bowtie2 -x /home/fgao/reference_genome/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome -U rep1_trim.fastq.gz -S rep1.sam -p 32
awk '{if($_!~"XS:i:") print }' rep1.sam > rep1.uniq.sam
samtools view -bS rep1.uniq.sam > rep1.uniq.bam
rm rep1.sam
rm rep1.uniq.sam
samtools sort rep1.uniq.bam -o rep1.sort.bam
samtools rmdup -s rep1.sort.bam rep1.rmdup.bam
samtools index rep1.rmdup.bam
