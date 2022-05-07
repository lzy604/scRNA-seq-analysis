## scRNA-seq-analysis(10x Genomics)

This pipeline performs the following tasks:
- Align reads, generate feature-barcode matrices ([Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorials))
- QC, analysis, and exploration([Seurat](https://satijalab.org/seurat/index.html)) 
- Running the pepline by MAESTRO


### System requirements
- Linux/Unix
- R 


### [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorials)
#### Installation
```bash
wget -O cellranger-6.1.2.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-6.1.2.tar.gz?Expires=1641366506&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci02LjEuMi50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2NDEzNjY1MDZ9fX1dfQ__&Signature=kaV8~ZabHhyDykUhbN~F78PDQfNZ64IamgsGc1nOSghFKPr0fbZ3WJk-2eWYh7IEt-KupenYP89W1zHi4lrxF~ZBbuP4NTaKEAa-G6ILJoX-VdyFnktkXFYDHgzEJ8ABq-NM6RWn20WD3a9BITNHTIWPtxjM-NaXAuR5uc5PuAEgjSDaQ2QBAQr~1q4aSM-~vJt~ia5e8acTz9RlM24EluLqfO59VCtAorP-5iJRwvLw9DjfrTlDtWfy3M2LSXp5OGmVJH1WUQReLK~0iZX2e8~vrHlAYpuxMa0Lgil6oHQ5s6vc~Dod3Aqpjb9sM~wuVo80zi4EqJ5nq0LU8SNbiQ__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
tar -zxvf cellranger-6.1.2.tar.gz
path=$(dirname cellranger-6.1.2)
export PATH=$path/cellranger-6.1.2:$PATH
cellranger
```
#### Running cellranger mkfastq(bcl2fastq2)
```bash
#rawdata download
wget https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-1.2.0.tar.gz
wget https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-simple-1.2.0.csv
tar -zxvf cellranger-tiny-bcl-1.2.0.tar.gz
path=$(dirname cellranger-tiny-bcl-simple-1.2.0.csv)
#mkfastq
cellranger mkfastq --id=project_name \
  --run=$path/cellranger-tiny-bcl-1.2.0 \
  --csv=$path/cellranger-tiny-bcl-simple-1.2.0.csv
```
fastqfiles are in the project_name/outs/fastq_path
#### Running cellranger count(aligns sequencing reads)
```bash
#reference transcriptome
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-3.0.0.tar.gz
tar -zxvf refdata-cellranger-GRCh38-3.0.0.tar.gz
refpath=$(dirname refdata-cellranger-GRCh38-3.0.0)
#fastq files download
wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar
tar -xvf pbmc_1k_v3_fastqs.tar
fastq_path=$(dirname pbmc_1k_v3_fastqs)
#count
cellranger count --id=project_ \
--fastqs= $fastq_path/pbmc_1k_v3_fastqs \
--sample=pbmc_1k_v3 \
--transcriptome= $refpath/refdata-cellranger-GRCh38-3.0.0
```
Outputs are in the pipestance directory in the outs folder.Load the filtered_feature_bc_matrix into third-party tools, such as Seurat.


### [Seurat](https://satijalab.org/seurat/index.html)
#### Installation
To install Seurat, R version 4.0 or greater is required. 
```R
# Enter commands in R (or R studio, if installed)
install.packages('Seurat')
library(Seurat)
```
#### Setup the Seurat Object
PBMC3k as the example
```bash
#download the data:
wget https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
tar -zxvf pbmc3k_filtered_gene_bc_matrices.tar.gz
```
Opean R
```

library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "../data/pbmc3k/filtered_gene_bc_matrices/hg19/")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
```
#### QC and selecting cells for further analysis
```R
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter cells that have unique feature counts over 2,500 or less than 200;have >5% mitochondrial counts
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
```

#### Identification of highly variable features (feature selection)
```
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
```

#### Scaling the data
```

``




## Running the pepline by MAESTRO
### Installation
```
#MAESTRO uses the Miniconda3 package management system to harmonize all of the software packages. 
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
conda config --add channels defaults
conda config --add channels liulab-dfci
conda config --add channels bioconda
conda config --add channels conda-forge
# To make the installation faster, we recommend using mamba
conda install mamba -c conda-forge
mamba create -n MAESTRO maestro=1.5.1 -c liulab-dfci
# Activate the environment
conda activate MAESTRO

# Download index files for GRCh38
mkdir MAESTRO/references
cd MAESTRO/references/
wget http://cistrome.org/~galib/MAESTRO/references/scRNA/Refdata_scRNA_MAESTRO_GRCh38_1.2.2.tar.gz
tar xvzf Refdata_scRNA_MAESTRO_GRCh38_1.2.2.tar.gz

# Build index for chromap. Only take a few minutes.
chromap -i -r Refdata_scATAC_MAESTRO_GRCh38_1.1.0/GRCh38_genome.fa -o GRCh38_chromap.index

```

### Example for 10x scRNA-seq analysis
```
# Before running MAESTRO, users need to activate the MAESTRO environment.
conda activate MAESTRO

# Please download the raw data from 10X genomics website.
wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar
tar xvf pbmc_1k_v3_fastqs.tar

# Configure the MAESTRO workflow
MAESTRO scrna-init --platform 10x-genomics --species GRCh38 --cores 32 --rseqc --directory pbmc_1k_v3_fastqs --count-cutoff 1000 --gene-cutoff 500 --cell-cutoff 10 --mapindex references/Refdata_scRNA_MAESTRO_GRCh38_1.2.2/GRCh38_STAR_2.7.6a --whitelist references/whitelist/3M-february-2018.txt --barcode-start 1 --barcode-length 16 --umi-start 17 --umi-length 12 --lisadir references/hg38_1000_2.0.h5 --signature human.immune.CIBERSORT

# Please download the raw data from 10X genomics website.
cd pbmc_1k_v3_fastqs/
wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar
tar xvf pbmc_1k_v3_fastqs.tar

# Configure the samples.json file
MAESTRO samples-init --assay_type scrna --platform 10x-genomics --data_typ fastq --data_dir data/pbmc_1k_v3_fastqs/

# Run snakemake pipeline
snakemake -np
nohup snakemake -j 16 > run.out &
```
    
