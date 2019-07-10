#!/bin/bash

set -e
set -o pipefail
set -o verbose

folder=`pwd`/tools/
cd $folder

#download STAR
wget https://github.com/alexdobin/STAR/archive/2.7.0f.tar.gz
tar xzf 2.7.0f.tar.gz
rm 2.7.0f.tar.gz

#download STAR-Fusion
wget https://github.com/STAR-Fusion/STAR-Fusion/releases/download/v1.6.0/STAR-Fusion-v1.6.0.FULL.tar.gz
tar xzf STAR-Fusion-v1.6.0.FULL.tar.gz
rm STAR-Fusion-v1.6.0.FULL.tar.gz

#download igv-reports, neccessary for STAR-fusion
wget https://github.com/igvteam/igv-reports/archive/v0.9.1.tar.gz
tar xzf v0.9.1.tar.gz
rm v0.9.1.tar.gz
cd igv-reports-0.9.1
python3 setup.py build

#download subread (featureCounts)
wget http://downloads.sourceforge.net/project/subread/subread-1.6.4/subread-1.6.4-Linux-x86_64.tar.gz
tar xzf subread-1.6.4-Linux-x86_64.tar.gz
rm subread-1.6.4-Linux-x86_64.tar.gz
mv subread-1.6.4-Linux-x86_64 subread-1.6.4

#download skewer (single-end adapter trimming)
git clone --depth 1 --branch 0.2.2 https://github.com/relipmoc/skewer.git skewer-0.2.2
cd skewer-0.2.2
make

#download salmon
wget https://github.com/COMBINE-lab/salmon/releases/download/v0.13.1/salmon-0.13.1_linux_x86_64.tar.gz
tar xzf salmon-0.13.1_linux_x86_64.tar.gz
rm salmon-0.13.1_linux_x86_64.tar.gz
mv salmon-latest_linux_x86_64 salmon-0.13.1
