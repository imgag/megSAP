#!/bin/bash
set -e
set -o pipefail
set -o verbose

root=`pwd`
folder=$root/tools/
cd $folder

#Ignore this - used for local installation
#folder=/mnt/storage2/megSAP/tools/

# #download hla-genotyper
cd $folder
wget https://github.com/axelgschwind/hla-genotyper/archive/refs/tags/2022_05.tar.gz
tar xzf 2022_05.tar.gz
rm 2022_05.tar.gz

#download varscan2
cd $folder
mkdir VarScan.v2.4.6
wget https://github.com/dkoboldt/varscan/raw/master/VarScan.v2.4.6.jar -O VarScan.v2.4.6/VarScan.v2.4.6.jar

#download strelka2
cd $folder
wget https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2
tar xjf strelka-2.9.10.centos6_x86_64.tar.bz2
rm strelka-2.9.10.centos6_x86_64.tar.bz2

#download mantis
cd $folder
wget https://github.com/OSU-SRLab/MANTIS/archive/v1.0.5.tar.gz
tar xzf v1.0.5.tar.gz
rm v1.0.5.tar.gz
cd MANTIS-1.0.5/tools/
make

#download umiVar2
cd $folder
git clone https://github.com/imgag/umiVar2.git umiVar2_2024_07
cd umiVar2_2024_07
git checkout 2024_07
cd ..
$folder/R-4.1.0/bin/R -f $root/install_deps_umiVar2.R

#install scarHRD
cd $folder
git clone https://github.com/imgag/scarHRD.git
$folder/R-4.1.0/bin/R -f $root/install_deps_scarhrd.R
