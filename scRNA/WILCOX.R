library(Seurat)
library(Matrix)
library(Matrix.utils)

data_cb<-Read10X(data.dir = "/home/fgao/Data_single_cell/Cai_Lab/GSE165371_cb_adult_mouse", gene.column=1)

#use genes expressed in 0.5% cells
sample_cb <- CreateSeuratObject(counts = data_cb, project = "CB")
meta<-read.csv("/home/fgao/Data_single_cell/Cai_Lab/GSE165371_cb_adult_mouse/cerebellum-atlas-analysis/data/cluster_metadata.csv", header=T)

# "Granule", "Bergmann", "MLI1", and "Purkinje" DEG analysis
sample_cb$cluster<-meta$cluster
Idents(sample_cb)<-meta$region

setwd("/home/fgao/Data_single_cell/Cai_Lab/GSE165371_cb_adult_mouse/")

library(Scillus)
deg1 <- find_diff_genes(dataset = sample_cb, clusters = c("AN1", "AN2", "COP", "CUL", "F", "I", "II", "III", "IX", "PF", "PRM", "SIM", "VI", "VII", "VIII", "X"), comparison = c("cluster", "Granule", "Bergmann"), logfc.threshold = 0, min.cells.group = 1, test.use = "wilcox")   
write.table(deg1, file="wilcox_Granule_Bergmann.txt", quote=F, sep="\t")

deg2 <- find_diff_genes(dataset = sample_cb, clusters = c("AN1", "AN2", "COP", "CUL", "F", "I", "II", "III", "IX", "PF", "PRM", "SIM", "VI", "VII", "VIII", "X"), comparison = c("cluster", "Granule", "MLI1"), logfc.threshold = 0, min.cells.group = 1, test.use = "wilcox")   
write.table(deg2, file="wilcox_Granule_MLI1.txt", quote=F, sep="\t")

deg3 <- find_diff_genes(dataset = sample_cb, clusters = c("AN1", "AN2", "COP", "CUL", "F", "I", "II", "III", "IX", "PF", "PRM", "SIM", "VI", "VII", "VIII", "X"), comparison = c("cluster", "Granule", "Purkinje"), logfc.threshold = 0, min.cells.group = 1, test.use = "wilcox")   
write.table(deg3, file="wilcox_Granule_Purkinje.txt", quote=F, sep="\t")

deg4 <- find_diff_genes(dataset = sample_cb, clusters = c("AN1", "AN2", "COP", "CUL", "F", "I", "II", "III", "IX", "PF", "PRM", "SIM", "VI", "VII", "VIII", "X"), comparison = c("cluster", "Bergmann", "MLI1"), logfc.threshold = 0, min.cells.group = 1, test.use = "wilcox")   
write.table(deg4, file="wilcox_Bergmann_MIL1.txt", quote=F, sep="\t")

deg5 <- find_diff_genes(dataset = sample_cb, clusters = c("AN1", "AN2", "COP", "CUL", "F", "I", "II", "III", "IX", "PF", "PRM", "SIM", "VI", "VII", "VIII", "X"), comparison = c("cluster", "Bergmann", "Purkinje"), logfc.threshold = 0, min.cells.group = 1, test.use = "wilcox")   
write.table(deg5, file="wilcox_Bergmann_Purkinje.txt", quote=F, sep="\t")

deg6 <- find_diff_genes(dataset = sample_cb, clusters = c("AN1", "AN2", "COP", "CUL", "F", "I", "II", "III", "IX", "PF", "PRM", "SIM", "VI", "VII", "VIII", "X"), comparison = c("cluster", "MLI1", "Purkinje"), logfc.threshold = 0, min.cells.group = 1, test.use = "wilcox")   
write.table(deg6, file="wilcox_MLI1_Purkinje.txt", quote=F, sep="\t")

