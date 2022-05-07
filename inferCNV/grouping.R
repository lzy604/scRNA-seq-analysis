setwd("/Users/liziyi/Documents/liu.data/sc-RNA,sc-ATAC/inferCNV/BRCA_GSE118389/")
##read the infercnv grouping
infercnv.group <- read.csv("./1k.10x.k2.c1/infercnv.observation_groupings.txt",sep = " ",row.names = 1)
infercnv.group2 <- data.frame(row.names = rownames(infercnv.group),infercnv=infercnv.group[,1])
group.1 <- rownames(infercnv.group2)[which(infercnv.group2$infercnv == 1)]
group.2 <- rownames(infercnv.group2)[which(infercnv.group2$infercnv == 2)]
##read the obj
infercnv_obj = readRDS('./1k.10x.k2.c1/run.final.infercnv_obj')
BRCA_GSE118389.cnv = as.data.frame(t(infercnv_obj@expr.data))
var.BRCA_GSE118389.cnv <- as.data.frame(apply(BRCA_GSE118389.cnv,1,function(x)var(x)))
var.BRCA_GSE118389.cnv.group1 <- var.BRCA_GSE118389.cnv[which(row.names(var.BRCA_GSE118389.cnv) %in% group.1),]
var.BRCA_GSE118389.cnv.group2 <- var.BRCA_GSE118389.cnv[which(row.names(var.BRCA_GSE118389.cnv) %in% group.2),]
if (mean(var.BRCA_GSE118389.cnv.group1)>mean(var.BRCA_GSE118389.cnv.group2)){
 tmp=data.frame(sample=group.1,anno=rep("cancer",length(group.1))) 
 tmp2=data.frame(sample=group.2,anno=rep("normal",length(group.2))) 
 infercnv.anno=rbind(tmp,tmp2)
}else{
  tmp=data.frame(sample=group.1,anno=rep("normal",length(group.1))) 
  tmp2=data.frame(sample=group.2,anno=rep("cancer",length(group.2))) 
  infercnv.anno=rbind(tmp,tmp2)
}
