setwd("/Users/liziyi/Documents/liu.data/sc-RNA,sc-ATAC/ENAH/HNSCC_GSE103322/")
library(Seurat)

SeuratObj <- readRDS("HNSCC_normal_SeuratObj.rds")
SeuratMarkers <- readRDS("HNSCC_normal_DiffMarkers.rds")
current.cluster.ids = as.integer(levels(SeuratMarkers$cluster))
SeuratObj@meta.data$assign.ident = SeuratObj@ident[rownames(SeuratObj@meta.data)]
manually.cluster.ids = c("0:Myofibroblasts","1:CD4T_Naive","2:CD8T","3:Fibroblasts","4:Fibroblasts","5:Treg_FOXP3","6:Endothelial","7:CD8T_Exhaust","8:Plasma","9:DC","10:Mast","11:Fibroblast","12:DC","13:Fibroblasts","14:T_MKI67","15:Endothelial","16:Myocyte")  
metadata <- SeuratObj@meta.data
SeuratObj@meta.data$assign.ident <- plyr::mapvalues(x = SeuratObj@meta.data$assign.ident, from = current.cluster.ids, to = manually.cluster.ids)
CAF.cell <- row.names(metadata)[grep("Fibroblast",SeuratObj@meta.data$assign.ident)]
CAF.cell.data <- SeuratObj@data[,colnames(SeuratObj@data) %in% CAF.cell]
cytokine <- read.csv("/Users/liziyi/Documents/liu.data/sc-RNA,sc-ATAC/ENAH/human.cytokine.txt",col.names = "cyto")[,1]
CAF.cell.cyto.data<- CAF.cell.data[row.names(CAF.cell.data) %in% cytokine,]
enah <- CAF.cell.data[row.names(CAF.cell.data) == "ENAH",]
HNSCC_GSE103322.cor <- as.data.frame(apply(CAF.cell.cyto.data,1,function(x)cor.test(enah,x)$p.value))
HNSCC_GSE103322.cor$cor <- apply(CAF.cell.cyto.data,1,function(x)cor.test(enah,x)$estimate)
colnames(HNSCC_GSE103322.cor)<- c("HNSCC_GSE103322.p.value","HNSCC_GSE103322.cor")

setwd("/Users/liziyi/Documents/liu.data/sc-RNA,sc-ATAC/ENAH/SKCM_GSE72056/")
SeuratObj <- readRDS("SKCM_normal_SeuratObj.rds")
SeuratMarkers <- readRDS("SKCM_normal_DiffMarkers.rds")
current.cluster.ids = as.integer(levels(SeuratMarkers$cluster))
SeuratObj@meta.data$assign.ident = SeuratObj@ident[rownames(SeuratObj@meta.data)]
manually.cluster.ids = c("0:CD4T_Memory","1:CD8T_Exhaust","2:BNaive","3:Myofibroblasts","4:CD8T_Exhaust","5:Monocyte","6:Fibroblasts","7:T_MKI67","8:Others","9:Monocyte","10:Treg_FOXP3","11:Endothelial")
SeuratObj@meta.data$assign.ident <- plyr::mapvalues(x = SeuratObj@meta.data$assign.ident, from = current.cluster.ids, to = manually.cluster.ids)
CAF.cell <- row.names(metadata)[grep("Fibroblast",SeuratObj@meta.data$assign.ident)]
CAF.cell.data <- SeuratObj@data[,colnames(SeuratObj@data) %in% CAF.cell]
cytokine <- read.csv("/Users/liziyi/Documents/liu.data/sc-RNA,sc-ATAC/ENAH/human.cytokine.txt",col.names = "cyto")[,1]
CAF.cell.cyto.data<- CAF.cell.data[row.names(CAF.cell.data) %in% cytokine,]
enah <- CAF.cell.data[row.names(CAF.cell.data) == "ENAH",]
SKCM_GSE72056.cor <- as.data.frame(apply(CAF.cell.cyto.data,1,function(x)cor.test(enah,x)$p.value))
SKCM_GSE72056.cor$cor <- apply(CAF.cell.cyto.data,1,function(x)cor.test(enah,x)$estimate)
colnames(SKCM_GSE72056.cor)<- c("SKCM_GSE72056.p.value","SKCM_GSE72056.cor")


setwd("/Users/liziyi/Documents/liu.data/sc-RNA,sc-ATAC/ENAH/NSCLC_EMTAB6149/")
SeuratObj <- readRDS("NSCLC_normal_SeuratObj.rds")
SeuratMarkers <- readRDS("NSCLC_normal_DiffMarkers.rds")
