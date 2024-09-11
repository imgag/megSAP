#!/bin/bash
set -e
set -o pipefail
set -o verbose

root=`pwd`
folder=`pwd`/tools/
python3_path=$folder/Python-3.10.9/

#Ignore this - used for local installation
#folder=/mnt/storage2/megSAP/tools/

#download minimap
cd $folder
curl -L https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 | tar -jxvf -

#download clair3
cd $folder
git clone https://github.com/HKU-BAL/Clair3.git Clair3-v1.0.10
cd Clair3-v1.0.10
git checkout "v1.0.10"
#download models
git clone https://github.com/nanoporetech/rerio.git
python3 rerio/download_model.py --clair3
mkdir -p models
mv rerio/clair3_models/* models
rm -rf rerio
cd models 
wget http://www.bio8.cs.hku.hk/clair3/clair3_models/r941_prom_sup_g5014.tar.gz
tar xzf r941_prom_sup_g5014.tar.gz
rm r941_prom_sup_g5014.tar.gz
wget http://www.bio8.cs.hku.hk/clair3/clair3_models/r1041_e82_400bps_sup_v420.tar.gz
tar xzf r1041_e82_400bps_sup_v420.tar.gz
rm r1041_e82_400bps_sup_v420.tar.gz
cd ..
#create temporary miniconda/clair env to build clair
wget https://repo.anaconda.com/miniconda/Miniconda3-py310_23.3.1-0-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p miniconda
source miniconda/bin/activate
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda create -p anaconda_clair3 -c bioconda clair3 python=3.9.0 -y
conda activate ./anaconda_clair3
#build clair3
make all
#cleanup
conda deactivate
conda deactivate # 2nd time to deactivate (base)
rm -rf anaconda_clair3/
rm miniconda.sh
rm -rf miniconda/
cd ..

#download pypy3
cd $folder
wget https://downloads.python.org/pypy/pypy3.9-v7.3.16-linux64.tar.bz2
tar xfvj pypy3.9-v7.3.16-linux64.tar.bz2
rm pypy3.9-v7.3.16-linux64.tar.bz2

#download parallel
cd $folder
wget https://ftpmirror.gnu.org/parallel/parallel-20240722.tar.bz2
tar xfvj parallel-20240722.tar.bz2
cd parallel-20240722
./configure --prefix=$folder/parallel-20240722
make
make install
cd ..
rm parallel-20240722.tar.bz2

#download longphase
cd $folder
mkdir -p longphase_v1.7.3
cd longphase_v1.7.3
wget https://github.com/twolinin/longphase/releases/download/v1.7.3/longphase_linux-x64.tar.xz
tar -xJf longphase_linux-x64.tar.xz
rm longphase_linux-x64.tar.xz
cd ..

#download straglr
cd $folder
wget https://github.com/bcgsc/straglr/releases/download/v1.5.1/straglr-1.5.1.tar.gz
tar xzf straglr-1.5.1.tar.gz 
rm straglr-1.5.1.tar.gz

#download StraglrOn
cd $folder
git clone https://github.com/imgag/StraglrOn.git StraglrOn_v0.2.2
cd StraglrOn_v0.2.2
git checkout "v0.2.2"
cd ..

#download Tandem Repeats Finder
cd $folder
mkdir -p TRF_4.09
cd TRF_4.09
wget https://github.com/Benson-Genomics-Lab/TRF/releases/download/v4.09.1/trf409.linux64
chmod +x trf409.linux64
ln -s trf409.linux64 trf

#download Blastn
cd $folder
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.16.0/ncbi-blast-2.16.0+-x64-linux.tar.gz
tar xzf ncbi-blast-2.16.0+-x64-linux.tar.gz 
rm ncbi-blast-2.16.0+-x64-linux.tar.gz 

#download modkit
cd $folder
wget https://github.com/nanoporetech/modkit/releases/download/v0.3.2/modkit_u16_x86_64.tar.gz
tar xzf modkit_u16_x86_64.tar.gz
mv dist_modkit_v0.3.2_828697b modkit_v0.3.2
chmod +x modkit_v0.3.2/modkit
rm modkit_u16_x86_64.tar.gz