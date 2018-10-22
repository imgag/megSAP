#!/bin/bash

set -e
set -o pipefail
set -o verbose

folder=`pwd`/tools/
cd $folder

#download STAR
git clone --depth 1 --branch 2.6.1a https://github.com/alexdobin/STAR.git STAR-2.6.1a
cd STAR-2.6.1a
git checkout 2.6.1a

#download STAR-Fusion
wget https://github.com/STAR-Fusion/STAR-Fusion/releases/download/STAR-Fusion-v1.4.0/STAR-Fusion-v1.4.0.FULL.tar.gz
tar xzf STAR-Fusion-v1.4.0.FULL.tar.gz
rm STAR-Fusion-v1.4.0.FULL.tar.gz

#download subread (featureCounts, subjunc)
wget http://downloads.sourceforge.net/project/subread/subread-1.6.2/subread-1.6.2-Linux-x86_64.tar.gz
tar xzf subread-1.6.2-Linux-x86_64.tar.gz
rm subread-1.6.2-Linux-x86_64.tar.gz
mv subread-1.6.2-Linux-x86_64 subread-1.6.2

#download skewer (single-end adapter trimming)
git clone --depth 1 --branch 0.2.2 https://github.com/relipmoc/skewer.git skewer-0.2.2
cd skewer-0.2.2
git checkout 0.2.2
make

#download salmon
wget https://github.com/COMBINE-lab/salmon/releases/download/v0.11.3/salmon-0.11.3-linux_x86_64.tar.gz
tar xzf salmon-0.11.3-linux_x86_64.tar.gz
rm salmon-0.11.3-linux_x86_64.tar.gz
mv salmon-0.11.3-linux_x86_64 salmon-0.11.3
