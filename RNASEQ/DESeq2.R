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
