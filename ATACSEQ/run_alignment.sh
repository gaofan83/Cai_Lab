chr_mouse="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX"
chr_human="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX"
chr_drosophila="chr2L chr2R chr3L chr3R chr4 chrX"

while read line;
do
bowtie2 -x /home/fgao/reference_genome/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome -1 ${line}_R1.fastq.gz -2 ${line}_R2.fastq.gz -S $line.sam -p 32
samtools view -bS $line.sam > $line.bam
rm $line.sam
samtools sort $line.bam -o $line.sort.bam
samtools rmdup -s $line.sort.bam $line.rmdup.bam
samtools idxstats $line.sort.bam > $line.mapping.txt
samtools idxstats $line.rmdup.bam > $line.mapping_rmdup.txt
done < sample_ID.txt

while read line;
do
samtools index $line.sort.bam
samtools index $line.rmdup.bam
samtools view -b $line.rmdup.bam $chr_mouse > $line.rmdup.nochrM.bam
samtools index $line.rmdup.nochrM.bam
samtools idxstats $line.sort.bam | awk '{s+=$3+$4} END {print "Total\t",s}' > $line.mapping.txt
samtools idxstats $line.sort.bam | awk '{s+=$3} END {print "Mapped\t",s}' >> $line.mapping.txt
samtools idxstats $line.sort.bam | awk '{if($1=="chrM") print "chrM\t",$3}' >> $line.mapping.txt
samtools idxstats $line.rmdup.bam | awk '{s+=$3} END {print "Mapped_rmdup\t",s}' >> $line.mapping.txt
samtools idxstats $line.rmdup.nochrM.bam | awk '{s+=$3} END {print "Mapped_rmdup_nochrM\t",s}' >> $line.mapping.txt
done < sample_ID.txt

 
export PATH=/home/fgao/anaconda3.5.2/bin:$PATH
export LD_LIBRARY_PATH=/home/fgao/anaconda3.5.2/lib:$LD_LIBRARY_PATH
export MANPATH=/home/fgao/anaconda3.5.2/share/man:$MANPATH
export PKG_CONFIG_PATH=/home/fgao/anaconda3.5.2/lib/pkgconfig:$PKG_CONFIG_PATH
while read line;
do
bamCoverage -b $line.rmdup.bam -o $line.bw --binSize 10 --normalizeUsing RPKM -p 16
done < sample_ID.txt
