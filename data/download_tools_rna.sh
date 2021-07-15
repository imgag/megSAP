#!/bin/bash

set -e
set -o pipefail
set -o verbose

folder=`pwd`/tools/
cd $folder

#download STAR
wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.9a.tar.gz
tar xzf 2.7.9a.tar.gz
rm 2.7.9a.tar.gz

#download STAR-Fusion
wget https://github.com/STAR-Fusion/STAR-Fusion/releases/download/v1.10.0/STAR-Fusion-v1.10.0.FULL.tar.gz
tar xzf STAR-Fusion-v1.10.0.FULL.tar.gz
rm STAR-Fusion-v1.10.0.FULL.tar.gz

#download igv-reports, neccessary for STAR-fusion
wget https://github.com/igvteam/igv-reports/archive/v0.9.1.tar.gz
tar xzf v0.9.1.tar.gz
rm v0.9.1.tar.gz
cd igv-reports-0.9.1 && python3 setup.py build && cd ..

#download subread (featureCounts)
wget http://downloads.sourceforge.net/project/subread/subread-2.0.2/subread-2.0.2-Linux-x86_64.tar.gz
tar xzf subread-2.0.2-Linux-x86_64.tar.gz
rm subread-2.0.2-Linux-x86_64.tar.gz
mv subread-2.0.2-Linux-x86_64 subread-2.0.2

#download skewer (single-end adapter trimming)
git clone --depth 1 --branch 0.2.2 https://github.com/relipmoc/skewer.git skewer-0.2.2
cd skewer-0.2.2 && make && cd ..
