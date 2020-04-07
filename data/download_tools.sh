#!/bin/bash

root=`pwd`
folder=$root/tools/

#download and build ngs-bits
cd $folder
git clone https://github.com/imgag/ngs-bits.git
cd ngs-bits
git checkout 2020_03 && git submodule update --recursive --init
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
git checkout v1.0.1 && git submodule update --recursive --init
make

#download ABRA2
cd $folder
mkdir abra2-2.22
cd abra2-2.22
wget https://github.com/mozack/abra2/releases/download/v2.22/abra2-2.22.jar -O abra2.jar

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

#download ClinCNV
cd $folder
git clone https://github.com/imgag/ClinCNV.git
cd ClinCNV
git fetch && git fetch --tags
git checkout 1.16.6
cd ..
mv ClinCNV ClinCNV-1.16.6

#download and build VEP
cd $root
chmod 755 download_tools_vep.sh
./download_tools_vep.sh

#download manta
cd $folder
wget https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2
tar xjf manta-1.6.0.centos6_x86_64.tar.bz2
rm manta-1.6.0.centos6_x86_64.tar.bz2
cd manta-1.6.0.centos6_x86_64
sed -i 's#referenceFasta = /illumina/development/Isis/Genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa##g' bin/configManta.py.ini

#download InterOp
cd $folder
wget https://github.com/Illumina/interop/releases/download/v1.0.25/InterOp-1.0.25-Linux-GNU-4.8.2.tar.gz
tar xzf InterOp-1.0.25-Linux-GNU-4.8.2.tar.gz
rm InterOp-1.0.25-Linux-GNU-4.8.2.tar.gz

#download Circos
cd $folder
wget http://circos.ca/distribution/circos-0.69-9.tgz
tar xzf circos-0.69-9.tgz
rm circos-0.69-9.tgz
circos-0.69-9/bin/circos -modules
