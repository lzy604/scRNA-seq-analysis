setwd("/Users/liziyi/Documents/liu.data/sc-RNA,sc-ATAC/inferCNV/BRCA_GSE118389/")
library("infercnv")
library("Seurat") 
#Running InferCNV
##simple 2-step protocol
###Creating an InferCNV object based on your three required inputs: the read count matrix, cell type annotations, and the gene ordering file:

#raw counts matrix
SeuratObj <- readRDS("BRCA_GSE118389_SeuratObj.rds")                                                                                          
counts_matrix = as.matrix(SeuratObj@data[,SeuratObj@cell.names])
dq.count.matrix <- (2^counts_matrix-1)*10
write.table(dq.count.matrix,file = "dq.10x.txt",sep = "\t",quote = FALSE)
#cell info
cell.info <- data.frame(row.names = colnames(dq.count.matrix))
cell.info$anno <- "cells"
write.table(cell.info,file = "cell.anno.txt",col.names = FALSE,sep = "\t",quote = FALSE)
#gene order file
# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix= "dq.10x.txt",
                                    annotations_file= "cell.anno.txt",
                                    delim="\t",
                                    gene_order_file= "gene_order_file.txt",
                                    ref_group_names=NULL)
# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="dq.k2.c1",  # dir is auto-created for storing outputs
                             cluster_by_groups=FALSE, 
                             k_obs_groups = 2,# cluster
                             denoise=T)





