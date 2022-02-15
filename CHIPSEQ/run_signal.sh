while read line;
do
./bigWigAverageOverBed ${line}.bw genome_bin.bed ${line}.tab
done < sample_ID.txt

./bigWigAverageOverBed ENCFF806NDV.bigWig genome_bin.bed ENCFF806NDV.tab
