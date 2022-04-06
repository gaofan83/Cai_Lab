# This is a pipeline to obtain complete MERVL-int sites across the genome.
# Repetitive element locations were downloaded from the RepeatMasker website (versions: mm10=4.0.5). 
# MERVL elements in particular were classified as “complete” if its internal part (MERVL-int) was flanked by two MT2_Mm elements facing in the same direction, separated by less than 7kb. 
# Original paper https://www.biorxiv.org/content/10.1101/523712v1.full.pdf

wget http://www.repeatmasker.org/genomes/mm10/RepeatMasker-rm405-db20140131/mm10.fa.out.gz
gunzip mm10.fa.out.gz
sed 's/[ ][ ]*/\t/g' mm10.fa.out > mm10.fa.out_reformat
awk '{if($10~"MERVL")print $5"\t"$6"\t"$7"\t"$10"\t1\t"$9}' mm10.fa.out_reformat > mm10_MERVL.bed
awk '{if($10~"MT2_Mm")print $5"\t"$6"\t"$7"\t"$10"\t1\t"$9}' mm10.fa.out_reformat > mm10_MT2_Mm.bed
sed -i 's/C/-/g' mm10_MERVL.bed
sed -i 's/C/-/g' mm10_MT2_Mm.bed
sortBed -i mm10_MERVL.bed > mm10_MERVL_sort.bed
sortBed -i mm10_MT2_Mm.bed > mm10_MT2_Mm_sort.bed
bedtools closest -a mm10_MERVL_sort.bed -b mm10_MT2_Mm_sort.bed -s -D ref -fu > mm10_MERVL_up_MT2_Mm.bed
bedtools closest -a mm10_MERVL_sort.bed -b mm10_MT2_Mm_sort.bed -s -D ref -fd > mm10_MERVL_dn_MT2_Mm.bed
paste mm10_MERVL_up_MT2_Mm.bed mm10_MERVL_dn_MT2_Mm.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$20"\t"$21"\t"$22"\t"$23"\t"$24"\t"$25"\t"$26}' > mm10_MERVL_up_dn_MT2_Mm.txt
awk '{if($13<0 && $20>0) print $0"\t"(($15+$16)/2-($8+$9)/2)}' mm10_MERVL_up_dn_MT2_Mm.txt > mm10_MERVL_up_dn_MT2_Mm_sel.txt
awk '{if($21<=7000) print}' mm10_MERVL_up_dn_MT2_Mm_sel.txt > mm10_MERVL_up_dn_MT2_Mm_7kb.txt

