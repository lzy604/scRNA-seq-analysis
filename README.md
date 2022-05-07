## scRNA-seq-analysis(10x Genomics)

This pipeline performs the following tasks:
- Align reads, generate feature-barcode matrices ([Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorials))
- QC, analysis, and exploration([Seurat](https://satijalab.org/seurat/index.html)) 
- Running the pepline by MAESTRO


### System requirements
- Linux/Unix
- Python
- R 


## Installation
We uses the Miniconda3 package management system to harmonize all of the software packages. 
Use the following commands to install Minicoda3：
``` bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
#### Create an isolated environment for RNA-seq
``` bash
conda create -n rna-seq
conda activate rna-seq
``` 

#### Install tools
Tools needed for this analysis are: R, samtools, FastQC, Trim Galore, STAR, RSeQC, stringtie, gffcompare, htseq-count. 
``` bash
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -c r r 
conda install -c bioconda samtools
conda install -c bioconda fastqc
conda install trim-galore
conda install STAR
conda install -c bioconda rseqc 
conda install -c bioconda htseq
conda install -c bioconda bioconductor-deseq2
conda install -c bioconda stringtie 
```

#### Genome files
Obtain a reference genome from Ensembl, iGenomes, NCBI or UCSC. In this example analysis we will use the mouse mm10 version of the genome from UCSC.
```bash
mkdir anno
cd anno
mkdir mm10
cd mm10
wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz

```

Generate genome indexes files for STAR mapping
```bash
tar zxvf Mus_musculus_UCSC_mm10.tar.gz
STAR --runThreadN 30 --runMode genomeGenerate --genomeDir star_index_mm10 --genomeFastaFiles /Sequence/WholeGenomeFasta/genome.fa --sjdbGTFfile /Annotation/Archives/archive-current/Genes/genes.gtf 
```
#### Download the public data

```
for ((i = 12;i<=15;i++)); #
do
fastq-dump --split-3 -O data/ SRR0020$i.sra.sra 
done
```

## Quality control on FastQ files 
FastQC aims to provide a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines. 

run FastQC interactively or using ht CLI, which offers the following options:
```bash
fastqc seqfile1 seqfile2 .. seqfileN
```

## Adapter Trim[OPTIONAL]
Use trim_glore to trim sequence adapter from the read FASTQ files.
```bash
trim_galore -q 20 --phred33 --stringency 3 --length 20 -e 0.1 \
            --paired $dir/cmp/01raw_data/$fq1 $dir/cmp/01raw_data/$fq2  \
            --gzip -o $input_data
```



## Running the pepline by MAESTRO
### installation
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
    
