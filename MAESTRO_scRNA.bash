
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

