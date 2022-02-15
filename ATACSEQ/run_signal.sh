while read line;
do
../../ENCODE/bigWigAverageOverBed ${line}.bw ../../ENCODE/genome_bin.bed ${line}.tab
done < sample_ID.txt
