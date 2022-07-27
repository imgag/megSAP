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

#download subread (featureCounts)
wget http://downloads.sourceforge.net/project/subread/subread-2.0.3/subread-2.0.3-Linux-x86_64.tar.gz
tar xzf subread-2.0.3-Linux-x86_64.tar.gz
rm subread-2.0.3-Linux-x86_64.tar.gz
mv subread-2.0.3-Linux-x86_64 subread-2.0.3

#download arriba
wget https://github.com/suhrig/arriba/releases/download/v2.2.1/arriba_v2.2.1.tar.gz
tar -xzf arriba_v2.2.1.tar.gz
rm arriba_v2.2.1.tar.gz
cd arriba_v2.2.1 && make && cd ..
conda create -c conda-forge -c bioconda -p arriba_v2.2.1/conda_env arriba=2.2.1