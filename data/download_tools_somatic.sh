#!/bin/bash
set -e
set -o pipefail
set -o verbose

root=`pwd`
folder=$root/tools/
cd $folder

#download hla-genotyper
wget https://github.com/axelgschwind/hla-genotyper/archive/refs/tags/2022_05.tar.gz
tar xzf 2022_05.tar.gz

#download varscan2
cd $folder
mkdir VarScan.v2.4.4
wget https://github.com/dkoboldt/varscan/raw/master/VarScan.v2.4.4.jar -O VarScan.v2.4.4/VarScan.v2.4.4.jar

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

#download umiVar2
cd $folder
git clone https://github.com/imgag/umiVar2.git umiVar2
cd umiVar2
# create virtual environment for umiVar2
# (requires installed python version from download_tools.sh)
$folder/Python3/bin/python3 -m venv venv 
source venv/bin/activate
pip3 install pysam==0.16.0.1
pip3 install numpy==1.19.5
pip3 install scipy==1.5.4
pip3 install networkx==2.5
deactivate
# install required R packages for umiVar2
# (requires installed R version 4.1.0 from download_tools.sh)
$folder/R-4.1.0/bin/R -f $root/install_deps_umiVar2.R