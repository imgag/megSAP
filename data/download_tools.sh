#!/bin/bash
set -e
set -o pipefail
set -o verbose

root=`pwd`
folder=$root/tools/
cpan_dir=$folder/perl_cpan/

#Ignore this - used for local installation
#folder=/mnt/storage2/megSAP/tools/
#cpan_dir=/mnt/storage2/megSAP/tools/perl_cpan/

python3_path=$folder/Python-3.10.9/

#download and build ngs-bits
cd $folder
git clone https://github.com/imgag/ngs-bits.git
cd ngs-bits
git checkout 2024_08 && git submodule update --recursive --init
make build_3rdparty
make build_tools_release

#download and build plain python3
cd $folder
mkdir -p $python3_path
cd $python3_path
wget https://www.python.org/ftp/python/3.10.9/Python-3.10.9.tgz
tar -zxvf Python-3.10.9.tgz
cd Python-3.10.9
./configure --prefix=$folder/Python-3.10.9 --enable-loadable-sqlite-extensions
make
make install

# create common python venv for megSAP
cd $folder
$folder/Python-3.10.9/bin/python3 -m venv Python-3.10.9_2024.08.21
source $folder/Python-3.10.9_2024.08.21/bin/activate
pip install --upgrade setuptools wheel
pip install -r $root/install_deps_python.txt --require-virtualenv
deactivate
cd ..

#download InterOp
cd $folder
wget https://github.com/Illumina/interop/releases/download/v1.2.4/interop-1.2.4-Linux-GNU.tar.gz
tar xzf interop-1.2.4-Linux-GNU.tar.gz
rm interop-1.2.4-Linux-GNU.tar.gz
