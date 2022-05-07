setwd("/project/dev/chenfei/project/SCRNAseq/Data_cancer/GBM_GSE84465/Result/infercnv/")
library("infercnv")
library(Seurat)
#Running InferCNV
##Creating an InferCNV object based on your three required inputs: the read count matrix, cell type annotations, and the gene ordering file:

###10x raw counts matrix
SeuratObj <- readRDS("GBM_GSE84465_SeuratObj.rds")                                                                                          
counts_matrix = as.matrix(SeuratObj@data[,SeuratObj@cell.names])
count.matrix <- (2^counts_matrix-1)*10
write.table(count.matrix,file = "10x.txt",sep = "\t",quote = FALSE)
##cell.anno
cell.info <-  data.frame(row.names = colnames(count.matrix))
cell.info$anno <-"cells"
write.table(cell.info,file = "cell.anno.txt",col.names = FALSE,sep = "\t",quote = FALSE)
##gene order file
# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix= "10x.txt",
                                    annotations_file= "cell.anno.txt",
                                    delim="\t",
                                    gene_order_file= "gene_order_file.txt",
                                    ref_group_names=NULL)
# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1,  # use 1 for smart-seq/Fluidigm C1, 0.1 for 10x-genomics/SNRS/MARS-Seq/Microwell/Drop-seq
                             out_dir="k2.c1",  # dir is auto-created for storing outputs
                             cluster_by_groups=FALSE,   # cluster
                             k_obs_groups = 2,
                             denoise=T)
##read the infercnv grouping
infercnv.group <- read.csv("./k2.c1/infercnv.observation_groupings.txt",sep = " ",row.names = 1)
infercnv.group <- data.frame(row.names = rownames(infercnv.group),infercnv=infercnv.group[,1])
group.1 <- rownames(infercnv.group)[which(infercnv.group$infercnv == 1)]
group.2 <- rownames(infercnv.group)[which(infercnv.group$infercnv == 2)]

##separate the normal and cancer cells
infercnv_obj = readRDS('./k2.c1/run.final.infercnv_obj')
infercnv = as.data.frame(t(infercnv_obj@expr.data))
var.infercnv <- as.data.frame(apply(infercnv,1,function(x)var(x)))
var.cnv.group1 <- var.infercnv[which(row.names(var.infercnv) %in% group.1),]
var.cnv.group2 <- var.infercnv[which(row.names(var.infercnv) %in% group.2),]
if (mean(var.cnv.group1)>mean(var.cnv.group2)){
  tmp=data.frame(sample=group.1,anno=rep("cancer",length(group.1))) 
  tmp2=data.frame(sample=group.2,anno=rep("normal",length(group.2))) 
  infercnv.anno=rbind(tmp,tmp2)
}else{
  tmp=data.frame(sample=group.1,anno=rep("normal",length(group.1))) 
  tmp2=data.frame(sample=group.2,anno=rep("cancer",length(group.2))) 
  infercnv.anno=rbind(tmp,tmp2)
}
write.table(infercnv.anno,"infercnv.anno.txt",col.names = F,row.names = FALSE,sep = "\t",quote = FALSE)

#normal cell
cnv.anno <- read.csv("infercnv.anno.txt",sep = "\t",row.names = 1,header = F)
colnames(cnv.anno)<- "infercnv"
normal.cell <- cnv.anno[(which(cnv.anno$infercnv == "normal")),]
#patient id
SeuratObj <- readRDS("GBM_GSE84465_SeuratObj.rds")
SeuratObj.info <- readLines("GSE84465_series_matrix.txt")
patient.info = SeuratObj.info[45]
patient.info = gsub('\"','',patient.info)
patient.info = gsub("patient id: ","",patient.info)
patient.info <- unlist(strsplit(patient.info,"\t"))[-1]

patient.info2 = SeuratObj.info[55]
patient.info2 = gsub('\"','',patient.info2)
patient.info2 <- unlist(strsplit(patient.info2,"\t"))[-1]

patient.id <- data.frame(patient.info,row.names = patient.info2)
row.names(patient.id)<- paste0("X",row.names(patient.id))

meta.data <- SeuratObj@meta.data
meta.data <- merge(meta.data,patient.id,0)
meta.data$orig.ident <-meta.data$patient.info
row.names(meta.data) <- meta.data$Row.names
meta.data <- meta.data[,-1]
meta.data <- merge(meta.data,cnv.anno,0)
row.names(meta.data) <- meta.data$Row.names
meta.data <- meta.data[,-1]
SeuratObj@meta.data <- meta.data
saveRDS(SeuratObj, file.path(paste0("/project/dev/chenfei/project/SCRNAseq/Data_cancer/GBM_GSE84465/",normal.SeuratObj@project.name, ".rds")))
pdf(file.path(paste0("/project/dev/chenfei/project/SCRNAseq/Data_cancer/GBM_GSE84465/Result/","tSNE_", SeuratObj@project.name, "_primary","_beforeCCA.pdf")),width=8, height=6)
TSNEPlot(object = SeuratObj, do.label = TRUE, pt.size = 0.5, group.by = "orig.ident")
dev.off()

##take out the normal cell
normal.SeuratObj <- SubsetData(SeuratObj, cells.use = rownames(SeuratObj@meta.data[which(SeuratObj@meta.data[,'infercnv']=='normal'),])); 
saveRDS(normal.SeuratObj, file.path(paste0("/project/dev/chenfei/project/SCRNAseq/Data_cancer/GBM_GSE84465/",normal.SeuratObj@project.name, "_infercnv_normal_beforeCCA_SeuratObj_.rds")))
pdf(file.path(paste0("/project/dev/chenfei/project/SCRNAseq/Data_cancer/GBM_GSE84465/","tSNE_", SeuratObj@project.name, "_infercnv_normal","beforeCCA.pdf")),width=8, height=6)
TSNEPlot(object = normal.SeuratObj, do.label = TRUE, pt.size = 0.5, group.by = "orig.ident")
dev.off()
