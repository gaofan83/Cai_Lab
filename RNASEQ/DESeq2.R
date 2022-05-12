library('DESeq2')
directory<-'./ft_seq/'
#use grep to search for the 'treated' part of filename to collect files
Files<-grep("*",list.files(directory),value=TRUE)

ES<-Files[c(1,2)]
NMuMG<-Files[c(5,6)]

sampleCondition1<-c(rep('ES',2), rep('NMuMG',2))

deg1<-c(ES, NMuMG)
sampleTable<-data.frame(sampleName=deg1, fileName=deg1, condition=sampleCondition1)
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~condition)
colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c('ES','NMuMG'))
dds<-DESeq(ddsHTSeq)
res<-results(dds)
write.table(res,file="DEG_ES_NMuMG_DESEQ2.txt",quote=FALSE,sep="\t")
resOrdered<-res[order(res$padj),]
res_sig = subset(resOrdered, padj<0.05)
write.table(res_sig,file="DEG_ES_NMuMG_DESEQ2_SIG.txt",quote=FALSE,sep="\t")


#generate normalized gene count matrix
sampleCondition<-c(rep('MuERVLZscan4',3), rep('NEG',3), rep('Zscan4',4))
sampleTable<-data.frame(sampleName=Files, fileName=Files, condition=sampleCondition)
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~condition)
dds<-DESeq(ddsHTSeq)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file= "deseq2_gene_counts_normalized.txt", sep="\t", quote=F, col.names=T)


#generate genome-bin associated DESeq2 table
genes<-read.table("features.tsv", header=F, row.names=1)
mm10_bins<-read.csv("mm10_25kb_bins_genes.csv", header=F)

deseq2_1<-read.table("DEG_MuERVLZscan4_NEG_DESEQ2.txt", header=T, row.names=1)
deseq2_1_gene<-merge(deseq2_1, genes, by=0)
mm10_bins_1<-merge(mm10_bins, deseq2_1_gene, by.x="V7", by.y="V2", all=TRUE)
mm10_bins_1_adj<-mm10_bins_1[!is.na(mm10_bins_1$Row.names), -c(8, 16, 17)]
colnames(mm10_bins_1_adj)<-c("Gene_Name", "Bin_ID", "Chr", "Start", "End", "Strand", "Score", "Ensembl_ID", "baseMean", "log2FoldChange ", "lfcSE", "stat", "pvalue", "padj")
write.table(mm10_bins_1_adj, file="mm10_25kb_bins_genes_MuERVLZscan4_NEG_DESEQ2.txt", quote=F, sep="\t", row.names=F)
