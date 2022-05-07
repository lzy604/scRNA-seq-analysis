library(MAESTRO)
library(Seurat)
library(ggplot2)
library(RColorBrewer)

# load data
pbmc.gene <- Read10X('Result/STAR/pbmc_1k_v3Solo.out/Gene/filtered')

# Perform clustering and differential gene analysis
pbmc.RNA.res <- RNARunSeurat(inputMat = pbmc.gene,
                             project = "pbmc_1k_v3",
                             orig.ident = NULL,
                             min.c = 10,#Cell coverage
                             min.g = 200,#Gene coverage
                             mito = TRUE,#Mitochondrial gene filter
                             mito.cutoff = 20,#Filter the cells with more than 20% mitochondrial reads
                             variable.genes = 2000,#Variance shrinkage and high-variance gene identification
                             organism = "GRCh38",
                             dims.use = 1:15,#Dimension reduction and determining significant components
                             cluster.res = 0.6,
                             only.pos = FALSE,
                             genes.test.use = "presto",
                             genes.cutoff = 1e-05,
                             genes.pct = 0.1,
                             genes.logfc = 0.25,
                             outdir = "./docs/assets/scrna/")

# Annotate cell types
data(human.immune.CIBERSORT)
pbmc.RNA.res$RNA <- RNAAnnotateCelltype(RNA = pbmc.RNA.res$RNA,
                                        gene = pbmc.RNA.res$genes,
                                        signatures = "human.immune.CIBERSORT",
                                        min.score = 0.05,
                                        outdir = "./docs/assets/scrna/")

DimPlot(object = pbmc.RNA.res$RNA, label = TRUE, pt.size = 0.15,
        group.by = "assign.ident", label.size = 3,
        cols = brewer.pal(8,"Set2")) +
  theme_linedraw() + NoLegend()

# Identify driver transcription regulators
pbmc.RNA.tfs <- RNAAnnotateTranscriptionFactor(RNA = pbmc.RNA.res$RNA,
                                               genes = pbmc.RNA.res$genes,
                                               project = pbmc.RNA.res$RNA@project.name,
                                               organism = "GRCh38",
                                               lisa.path = "/Users/galib/miniconda3/envs/lisa_env/bin/.venvs/lisa_env/bin/lisa",
                                               top.tf = 10,
                                               outdir = "./docs/assets/scrna/")

# Visualize driver transcription factors for each cluster
tfs <- sapply(pbmc.RNA.monocyte.tfs[[1]], function(x){
  return(unlist(strsplit(x, split = " | ", fixed = TRUE))[1])})
VisualizeTFenrichment(TFs = tfs,
                      cluster.1 = 3,
                      type = "RNA",
                      SeuratObj = pbmc.RNA.res$RNA,
                      LISA.table = "pbmc_1k_v3_Monocyte_lisa.txt",
                      visual.totalnumber = 100,
                      name = "pbmc_1k_v3_Monocyte_filtered")

# Visualize TF/genes expression level using Vlnplot and Umap.
VisualizeVlnplot(genes = c("TF1","TF2"),
                 type = "RNA",
                 SeuratObj = pbmc.RNA.res$RNA,
                 ncol = 2,
                 width = 8.5,
                 height = 4,
                 name = "pbmc_1k_v3_Monocyte")
VisualizeUmap(genes = c("TF1","TF2"),
              type = "RNA",
              SeuratObj = pbmc.RNA.res$RNA,
              ncol = 2,
              width = 8,
              height = 3,
              name = "pbmc_1k_v3_Monocyte")
# Save the object for future analysis
saveRDS(pbmc.RNA.res, "pbmc_1k_v3_8k_res.rds")


