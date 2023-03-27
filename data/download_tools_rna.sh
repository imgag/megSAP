#!/bin/bash
set -e
set -o pipefail
set -o verbose

root=`pwd`
folder=`pwd`/tools/

#Ignore this - used for local installation
#folder=/mnt/storage2/megSAP/tools/

#download STAR
cd $folder
wget https://github.com/alexdobin/STAR/archive/2.7.10b.tar.gz
tar xzf 2.7.10b.tar.gz
rm 2.7.10b.tar.gz

#download subread (featureCounts)
cd $folder
wget https://downloads.sourceforge.net/project/subread/subread-2.0.4/subread-2.0.4-Linux-x86_64.tar.gz
tar xzf subread-2.0.4-Linux-x86_64.tar.gz
rm subread-2.0.4-Linux-x86_64.tar.gz
mv subread-2.0.4-Linux-x86_64 subread-2.0.4

#download arriba
cd $folder
wget https://github.com/suhrig/arriba/releases/download/v2.4.0/arriba_v2.4.0.tar.gz
tar -xzf arriba_v2.4.0.tar.gz
rm arriba_v2.4.0.tar.gz
cd arriba_v2.4.0 && make && cd ..
#install R dependencies:
$folder/R-4.1.0/bin/R -f $root/install_deps_arriba.R

#download kraken2
cd $folder
wget https://github.com/DerrickWood/kraken2/archive/refs/tags/v2.1.2.tar.gz
tar -xzf v2.1.2.tar.gz
cd kraken2-2.1.2
./install_kraken2.sh bin
cd ..
rm v2.1.2.tar.gz

