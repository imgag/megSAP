#!/bin/bash
set -e
set -o pipefail
set -o verbose

root=`pwd`
#adapt if you want to store elsewhere
folder=$root/tools/singularity_container

#Ignore this - used for local installation
#folder=/mnt/storage2/megSAP/tools/singularity_container

# make sure container folder exists
mkdir -p $folder
cd $folder

#clair3-trio
singularity pull docker://hkubal/clair3-trio:v0.7
