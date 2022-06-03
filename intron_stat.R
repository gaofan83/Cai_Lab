install.packages("remotes")
remotes::install_github("asrinivasan-oa/gread")

library(gread)
library(Repitools)

gtf_file <- file.path("/home/fgao/refdata-cellranger-mm10-3.0.0/genes/genes.gtf")
gtf <- read_format(gtf_file)
ans <- construct_introns(gtf, update=TRUE)[] #add intron information
introns <- construct_introns(gtf, update=FALSE) #only intron information
mat_introns <- annoGR2DF(introns)
out_introns <- mat_introns[,c("width", "gene_name", "transcript_name")]
out_introns$count <- rep(1, nrow(out_introns))
length_introns <- aggregate(width ~ transcript_name, data=out_introns, FUN=sum) 
count_introns <- aggregate(count ~ transcript_name, data=out_introns, FUN=sum) 
id_introns <- unique(out_introns[, c("gene_name", "transcript_name")])

length_introns_genes <- merge(id_introns, length_introns, by="transcript_name")
count_introns_genes <- merge(id_introns, count_introns, by="transcript_name")

length_intron_genes_stat <- aggregate(width ~ gene_name, data=length_introns_genes, FUN=mean) 
count_intron_genes_stat <- aggregate(count ~ gene_name, data=count_introns_genes, FUN=mean) 

write.table(length_intron_genes_stat, file="gene_intron_length_mean.txt", quote=F, sep="\t", row.names=F)
write.table(count_intron_genes_stat, file="gene_intron_count_mean.txt", quote=F, sep="\t", row.names=F)
