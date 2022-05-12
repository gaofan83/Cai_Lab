awk -F '"' '{if ($2!="") print }' 2020-01-30-mouse-whole-genome-25kb-blocks.csv | sed 's/ /\n/g' | \
  awk -F '"'  '{if($1~"^chr") {pos=$1; print} else print pos$1}' | sed 's/"//g' | sed 's/\n,/\n/g' > \
  mm10_25kb_bins_genes.csv
  
