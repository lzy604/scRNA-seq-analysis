
# MAESTRO uses the Miniconda3 package management system to harmonize all of the software packages. 
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

# Download lisa2 data files for GRCh38
cd ../
mkdir annotation
cd annotation/
wget http://cistrome.org/~alynch/data/lisa_data/hg38_1000_2.0.h5

# Download giggle annotation files for GRCh38
wget http://cistrome.org/~galib/MAESTRO/references/giggle.all.tar.gz
tar xvzf giggle.all.tar.gz
