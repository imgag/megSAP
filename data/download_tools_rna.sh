#!/bin/bash

set -e
set -o pipefail
set -o verbose

folder=`pwd`/tools/
cd $folder

#download STAR
wget https://github.com/alexdobin/STAR/archive/2.7.10a.tar.gz
tar xzf 2.7.10a.tar.gz
rm 2.7.10a.tar.gz

#download STAR-Fusion
wget https://github.com/STAR-Fusion/STAR-Fusion/releases/download/v1.9.1/STAR-Fusion-v1.9.1.FULL.tar.gz
tar xzf STAR-Fusion-v1.9.1.FULL.tar.gz
rm STAR-Fusion-v1.9.1.FULL.tar.gz
sed -i '1s/python/python3/' STAR-Fusion-v1.9.1/FusionInspector/FusionInspector

#download samtools 1.7, this version is neccessary for STAR-fusion
wget https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2
tar xjf samtools-1.7.tar.bz2
rm samtools-1.7.tar.bz2
cd samtools-1.7
make

#download subread (featureCounts)
wget http://downloads.sourceforge.net/project/subread/subread-2.0.3/subread-2.0.3-Linux-x86_64.tar.gz
tar xzf subread-2.0.3-Linux-x86_64.tar.gz
rm subread-2.0.3-Linux-x86_64.tar.gz
mv subread-2.0.3-Linux-x86_64 subread-2.0.3

#download skewer (single-end adapter trimming)
git clone --depth 1 --branch 0.2.2 https://github.com/relipmoc/skewer.git skewer-0.2.2
cd skewer-0.2.2 && make && cd ..

#download arriba
wget https://github.com/suhrig/arriba/releases/download/v2.2.1/arriba_v2.2.1.tar.gz
tar -xzf arriba_v2.2.1.tar.gz
rm arriba_v2.2.1.tar.gz
cd arriba_v2.2.1 && make && cd ..
conda create -c conda-forge -c bioconda -p arriba_v2.2.1/conda_env arriba=2.2.1