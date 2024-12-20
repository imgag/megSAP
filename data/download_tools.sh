#!/bin/bash
set -e
set -o pipefail
set -o verbose

root=`pwd`
folder=$root/tools/

#Ignore this - used for local installation
#folder=/mnt/storage2/megSAP/tools/

#download InterOp
cd $folder
wget https://github.com/Illumina/interop/releases/download/v1.2.4/interop-1.2.4-Linux-GNU.tar.gz
tar xzf interop-1.2.4-Linux-GNU.tar.gz
rm interop-1.2.4-Linux-GNU.tar.gz
