#!/bin/bash
mkdir -p `pwd`/tools
folder=`pwd`/tools/

cd $folder

#download RepeatMasker
wget http://www.repeatmasker.org/RepeatMasker-open-4-0-6.tar.gz
tar xzf RepeatMasker-open-4-0-6.tar.gz
rm -rf RepeatMasker-open-4-0-6.tar.gz

#download and build ngs-bits
cd $folder
git clone https://github.com/imgag/ngs-bits.git --depth 1
cd ngs-bits
git checkout 2018_06 && git submodule update --recursive --init
make build_3rdparty
make build_tools_release

#download and build samtools
cd $folder
wget http://downloads.sourceforge.net/project/samtools/samtools/1.6/samtools-1.6.tar.bz2
tar xjf samtools-1.6.tar.bz2
rm samtools-1.6.tar.bz2
cd samtools-1.6
make

#download and build freebayes
cd $folder
git clone https://github.com/ekg/freebayes.git --depth 1
cd freebayes
git checkout v1.1.0 && git submodule update --recursive --init
make

#download and build vcflib
cd $folder
git clone https://github.com/vcflib/vcflib.git --depth 1
cd vcflib
git checkout v1.0.0-rc1 && git submodule update --recursive --init
make

#download ABRA2
cd $folder
mkdir abra2-2.15
cd abra2-2.15
wget https://github.com/mozack/abra2/releases/download/v2.15/abra2-2.15.jar -O abra2.jar

#download and build samblaster
cd $folder
git clone https://github.com/GregoryFaust/samblaster.git --depth 1
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

#download and build VEP
chmod 755 download_tools_vep.sh
./download_tools_vep.sh

#download and build ClinCnv
echo "TODO - including R"
