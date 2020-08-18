#!/bin/bash

set -e
set -o pipefail
set -o verbose

folder=`pwd`/tools/
cd $folder

#download STAR
wget https://github.com/alexdobin/STAR/archive/2.7.3a.tar.gz
tar xzf 2.7.3a.tar.gz
rm 2.7.3a.tar.gz

#download STAR-Fusion
wget https://github.com/STAR-Fusion/STAR-Fusion/releases/download/v1.9.1/STAR-Fusion-v1.9.1.FULL.tar.gz
tar xzf STAR-Fusion-v1.9.1.FULL.tar.gz
rm STAR-Fusion-v1.9.1.FULL.tar.gz
sed -i '1s/python/python3/' STAR-Fusion-v1.9.1/FusionInspector/FusionInspector

#download igv-reports, neccessary for STAR-fusion
wget https://github.com/igvteam/igv-reports/archive/v0.9.1.tar.gz
tar xzf v0.9.1.tar.gz
rm v0.9.1.tar.gz
cd igv-reports-0.9.1 && python3 setup.py build && cd ..

#download subread (featureCounts)
wget http://downloads.sourceforge.net/project/subread/subread-2.0.0/subread-2.0.0-Linux-x86_64.tar.gz
tar xzf subread-2.0.0-Linux-x86_64.tar.gz
rm subread-2.0.0-Linux-x86_64.tar.gz
mv subread-2.0.0-Linux-x86_64 subread-2.0.0

#download skewer (single-end adapter trimming)
git clone --depth 1 --branch 0.2.2 https://github.com/relipmoc/skewer.git skewer-0.2.2
cd skewer-0.2.2 && make && cd ..
