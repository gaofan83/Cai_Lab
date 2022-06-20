while read line;
do
./bigWigAverageOverBed ${line}.bw genome_bin.bed ${line}.tab
done < sample_ID.txt

./bigWigAverageOverBed ENCFF806NDV.bigWig genome_bin.bed ENCFF806NDV.tab

touch score_all.tab
while read line;
do
./bigWigAverageOverBed ${line}.bw genome_bin.bed ${line}.tab
paste score_all.tab ${line}.tab > score_new.tab
mv score_new.tab score_all.tab
done < sample_ID.txt

echo "CHR_ID,H3K27me3_rep1,H3K27me3_rep2,H3K27ac_rep1,H3K27ac_rep2,H3K4me1_rep1,H3K4me1_rep2,H3K4me3_rep1,H3K4me3_rep2,IgG_rep1,IgG_rep2" > WT_scores.txt
awk '{print $1","$6","$12","$18","$24","$30","$36","$42","$48","$54","$60}' score_all.tab >> WT_scores.txt

echo "POS,H3K27me3_rep1,H3K27me3_rep2,H3K27ac_rep1,H3K27ac_rep2,H3K4me1_rep1,H3K4me1_rep2,H3K4me3_rep1,H3K4me3_rep2,IgG_rep1,IgG_rep2" > GFP_scores.txt
awk '{print $1","$66","$72","$78","$84","$90","$96","$102","$108","$114","$120}' score_all.tab >> GFP_scores.txt
