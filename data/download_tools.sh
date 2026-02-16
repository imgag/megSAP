#!/bin/bash
set -e
set -o pipefail
set -o verbose

root=`dirname $(pwd)`

# Get data_folder path
SETTINGS_FILE=$root/settings.ini
if [ ! -f "$SETTINGS_FILE" ]; then
    SETTINGS_FILE="$root/settings.ini.default"
fi
DATA_FOLDER=$(grep -E "^data_folder" "$SETTINGS_FILE" | awk -F ' = ' '{print $2}' | sed "s|\[path\]|${root}|")

folder=$DATA_FOLDER/tools/

# Ensure the tools folder exists
mkdir -p $folder

#Ignore this - used for local installation
#folder=/mnt/storage2/megSAP/tools/

#download InterOp
cd $folder
wget https://github.com/Illumina/interop/releases/download/v1.2.4/interop-1.2.4-Linux-GNU.tar.gz
tar xzf interop-1.2.4-Linux-GNU.tar.gz
rm interop-1.2.4-Linux-GNU.tar.gz

#download dorado
cd $folder
wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-1.3.1-linux-x64.tar.gz
tar xzf dorado-1.3.1-linux-x64.tar.gz
rm dorado-1.3.1-linux-x64.tar.gz