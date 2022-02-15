export PATH=/home/fgao/anaconda3.5.2/bin:$PATH
export LD_LIBRARY_PATH=/home/fgao/anaconda3.5.2/lib:$LD_LIBRARY_PATH
export MANPATH=/home/fgao/anaconda3.5.2/share/man:$MANPATH
export PKG_CONFIG_PATH=/home/fgao/anaconda3.5.2/lib/pkgconfig:$PKG_CONFIG_PATH

#strand-specific signals
bamCoverage -b rep1.rmdup.bam -o rep1_fwd.bw --binSize 10 --normalizeUsing RPKM -p 32 --filterRNAstrand reverse
bamCoverage -b rep1.rmdup.bam -o rep1_rev.bw --binSize 10 --normalizeUsing RPKM -p 32 --filterRNAstrand forward

bamCoverage -b rep2.rmdup.bam -o rep2_fwd.bw --binSize 10 --normalizeUsing RPKM -p 32 --filterRNAstrand reverse
bamCoverage -b rep2.rmdup.bam -o rep2_rev.bw --binSize 10 --normalizeUsing RPKM -p 32 --filterRNAstrand forward

awk '{if($3=="gene" && $7=="+") print $1"\t"$4"\t"$5"\t"$14}' ~/gencode/gencode.vM25.annotation.gtf > mm10_gene_body_fwd.bed
sed -i 's/"//g' mm10_gene_body_fwd.bed
sed -i 's/;//g' mm10_gene_body_fwd.bed

awk '{if($3=="gene" && $7=="-") print $1"\t"$4"\t"$5"\t"$14}' ~/gencode/gencode.vM25.annotation.gtf > mm10_gene_body_rev.bed
sed -i 's/"//g' mm10_gene_body_rev.bed
sed -i 's/;//g' mm10_gene_body_rev.bed

computeMatrix scale-regions -S rep*fwd.bw -R mm10_gene_body_fwd.bed -a 0 -b 0 -m 10 --binSize 10 \
   -o matrix_scaled_fwd.gz --outFileNameMatrix matrix_scaled_fwd.tab --outFileSortedRegions mm10_gene_body_fwd_scores.bed

zcat matrix_scaled_fwd.gz | awk '{if($1~"chr") print}' > groseq_intensity_fwd.bed
paste groseq_intensity_fwd.bed mm10_gene_body_fwd.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$8"\t"$12}' > groseq_intensity_fwd.txt

computeMatrix scale-regions -S rep*rev.bw -R mm10_gene_body_rev.bed -a 0 -b 0 -m 10 --binSize 10 \
   -o matrix_scaled_rev.gz --outFileNameMatrix matrix_scaled_rev.tab --outFileSortedRegions mm10_gene_body_rev_scores.bed

zcat matrix_scaled_rev.gz | awk '{if($1~"chr") print}' > groseq_intensity_rev.bed
paste groseq_intensity_rev.bed mm10_gene_body_rev.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$8"\t"$12}' > groseq_intensity_rev.txt


awk '{if($3=="gene" && $7=="+") print $1"\t"$4-1"\t"$4"\t"$14;}' ~/gencode/gencode.vM25.annotation.gtf > mm10_tss_fwd.bed
awk '{if($3=="gene" && $7=="-") print $1"\t"$5"\t"$5+1"\t"$14;}' ~/gencode/gencode.vM25.annotation.gtf > mm10_tss_rev.bed

sed -i 's/"//g' mm10_tss_fwd.bed
sed -i 's/;//g' mm10_tss_fwd.bed
sed -i 's/"//g' mm10_tss_rev.bed
sed -i 's/;//g' mm10_tss_rev.bed

computeMatrix reference-point -S rep*fwd.bw -R mm10_tss_fwd.bed -a 500 -b 500 --binSize 500 -p 32 \
   -o matrix_tss_fwd_scaled.gz --outFileNameMatrix matrix_tss_fwd_scaled.tab --outFileSortedRegions mm10_tss_fwd_scores.bed

computeMatrix reference-point -S rep*rev.bw -R mm10_tss_rev.bed -a 500 -b 500 --binSize 500 -p 32 \
   -o matrix_tss_rev_scaled.gz --outFileNameMatrix matrix_tss_rev_scaled.tab --outFileSortedRegions mm10_tss_rev_scores.bed

zcat matrix_tss_fwd_scaled.gz | awk '{if($1~"chr") print}' > groseq_intensity_tss_fwd.bed
paste groseq_intensity_tss_fwd.bed mm10_tss_fwd.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"($7+$8)"\t"($9+$10)"\t"$14}' > groseq_intensity_tss_fwd.txt

zcat matrix_tss_rev_scaled.gz | awk '{if($1~"chr") print}' > groseq_intensity_tss_rev.bed
paste groseq_intensity_tss_rev.bed mm10_tss_rev.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"($7+$8)"\t"($9+$10)"\t"$14}' > groseq_intensity_tss_rev.txt
