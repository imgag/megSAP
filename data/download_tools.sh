#!/bin/bash

root=`pwd`
folder=$root/tools/

#download and build ngs-bits
cd $folder
git clone https://github.com/imgag/ngs-bits.git
cd ngs-bits
git checkout 2018_11 && git submodule update --recursive --init
make build_3rdparty
make build_tools_release

#download and build samtools
cd $folder
wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
tar xjf samtools-1.9.tar.bz2
rm samtools-1.9.tar.bz2
cd samtools-1.9
make

#download and build freebayes
cd $folder
git clone https://github.com/ekg/freebayes.git
cd freebayes
git checkout v1.2.0 && git submodule update --recursive --init
make

#download and build vcflib
cd $folder
git clone https://github.com/vcflib/vcflib.git
cd vcflib
git checkout v1.0.0-rc2 && git submodule update --recursive --init
make

#download ABRA2
cd $folder
mkdir abra2-2.18
cd abra2-2.18
wget https://github.com/mozack/abra2/releases/download/v2.18/abra2-2.18.jar -O abra2.jar

#download and build samblaster
cd $folder
git clone https://github.com/GregoryFaust/samblaster.git
cd samblaster
git checkout v.0.1.24
make

#download and build BWA
cd $folder
wget http://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.17.tar.bz2
tar xjf bwa-0.7.17.tar.bz2
rm bwa-0.7.17.tar.bz2
cd bwa-0.7.17
make

#download ClinCnv
cd $folder
git clone https://github.com/imgag/ClinCNV.git
cd ClinCNV
git checkout 1.07
cd ..
mv ClinCNV ClinCNV-1.07

#download and build VEP
cd $root
chmod 755 download_tools_vep.sh
./download_tools_vep.sh

