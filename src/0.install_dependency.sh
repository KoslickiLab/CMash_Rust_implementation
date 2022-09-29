#!/bin/bash

### Build conda reproducible ENVs
date

# active conda inside script
temp=$(which conda)
conda_path=$(echo ${temp%/*bin/conda})
if [ -f ${conda_path}/etc/profile.d/conda.sh ]; then
        . ${conda_path}/etc/profile.d/conda.sh
else
        echo "ERROR: conda path can't be corrected identified!!!"
        exit 1
fi
unset conda_path

# this is a local conda env (not in this repo)
conda create -y -n metagenome python=3.7
conda activate metagenome 
conda install -y -c anaconda seaborn
conda install -y -c bioconda kmc
conda install -y -c conda-forge -c bioconda sourmash
conda install -y -c bioconda sra-tools
conda install -y -c bioconda minimap2
conda install -y -c bioconda fastqc
conda install -y -c conda-forge marisa-trie
conda install -y -c conda-forge hydra
conda install -y -c conda-forge biopython
conda install -y -c bioconda khmer
conda install -y -c anaconda h5py
conda deactivate


