#!/bin/bash

root=`pwd`
folder=$root/tools/

#download and build ngs-bits
cd $folder
git clone https://github.com/imgag/ngs-bits.git
cd ngs-bits
git checkout 2022_04 && git submodule update --recursive --init
make build_3rdparty
make build_tools_release

#download and build samtools
cd $folder
wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2
tar xjf samtools-1.11.tar.bz2
rm samtools-1.11.tar.bz2
cd samtools-1.11
make

#download and build freebayes
cd $folder
git clone https://github.com/ekg/freebayes.git freebayes-1.3.3
cd freebayes-1.3.3
git checkout v1.3.3 && git submodule update --recursive --init
meson build/ --buildtype release
cd build
ninja
ninja test

#download and build vcflib
cd $folder
git clone https://github.com/vcflib/vcflib.git vcflib-1.0.2
cd vcflib-1.0.2
git checkout v1.0.2 && git submodule update --recursive --init
mkdir -p build && cd build
cmake ..
cmake --build .
cmake --install .

#download ABRA2
cd $folder
mkdir abra2-2.23
cd abra2-2.23
wget https://github.com/mozack/abra2/releases/download/v2.23/abra2-2.23.jar -O abra2.jar

#download and build samblaster
cd $folder
git clone https://github.com/GregoryFaust/samblaster.git
cd samblaster
git checkout v.0.1.26
make

#download and build BWA
cd $folder
wget http://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.17.tar.bz2
tar xjf bwa-0.7.17.tar.bz2
rm bwa-0.7.17.tar.bz2
cd bwa-0.7.17
make

#download bwa-mem2
cd $folder
wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.1/bwa-mem2-2.1_x64-linux.tar.bz2
tar xjf bwa-mem2-2.1_x64-linux.tar.bz2
rm bwa-mem2-2.1_x64-linux.tar.bz2

#download and build R (for ClinCnv and UmiVar2)
cd $folder
wget https://cran.r-project.org/src/base/R-4/R-4.1.0.tar.gz
tar -xvzf R-4.1.0.tar.gz
mv R-4.1.0 R-4.1.0-src
cd R-4.1.0-src
./configure --with-pcre1 --prefix $folder/R-4.1.0
make all install
cd ..
rm -rf R-4.1.0.tar.gz R-4.1.0-src

#download ClinCNV
cd $folder
git clone https://github.com/imgag/ClinCNV.git ClinCNV-1.17.2
mv ClinCNV-1.17.2/cytobandsHG38.txt ClinCNV-1.17.2/cytobands.txt
# install required R packages for ClinCNV
$folder/R-4.1.0/bin/R -f $root/install_deps_clincnv.R

#download and build VEP
cd $root
chmod 755 download_tools_vep.sh
./download_tools_vep.sh

#download manta
cd $folder
wget https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2
tar xjf manta-1.6.0.centos6_x86_64.tar.bz2
rm manta-1.6.0.centos6_x86_64.tar.bz2
cd manta-1.6.0.centos6_x86_64
sed -i 's#referenceFasta = /illumina/development/Isis/Genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa##g' bin/configManta.py.ini

#download InterOp
cd $folder
wget https://github.com/Illumina/interop/releases/download/v1.0.25/InterOp-1.0.25-Linux-GNU-4.8.2.tar.gz
tar xzf InterOp-1.0.25-Linux-GNU-4.8.2.tar.gz
rm InterOp-1.0.25-Linux-GNU-4.8.2.tar.gz

#download Circos
cd $folder
wget http://circos.ca/distribution/circos-0.69-9.tgz
tar xzf circos-0.69-9.tgz
rm circos-0.69-9.tgz

# install required Perl modules for Circos in a subfolder
circos_cpan_dir=$folder/circos-0.69-9/cpan/
mkdir -p $circos_cpan_dir
cpanm -l $circos_cpan_dir -L $circos_cpan_dir Carp
cpanm -l $circos_cpan_dir -L $circos_cpan_dir Clone
cpanm -l $circos_cpan_dir -L $circos_cpan_dir Config::General
cpanm -l $circos_cpan_dir -L $circos_cpan_dir Cwd
cpanm -l $circos_cpan_dir -L $circos_cpan_dir Data::Dumper
cpanm -l $circos_cpan_dir -L $circos_cpan_dir Digest::MD5
cpanm -l $circos_cpan_dir -L $circos_cpan_dir File::Basename
cpanm -l $circos_cpan_dir -L $circos_cpan_dir File::Spec::Functions
cpanm -l $circos_cpan_dir -L $circos_cpan_dir File::Temp
cpanm -l $circos_cpan_dir -L $circos_cpan_dir FindBin
cpanm -l $circos_cpan_dir -L $circos_cpan_dir Font::TTF::Font
cpanm -l $circos_cpan_dir -L $circos_cpan_dir GD
cpanm -l $circos_cpan_dir -L $circos_cpan_dir GD::Polyline
cpanm -l $circos_cpan_dir -L $circos_cpan_dir Getopt::Long
cpanm -l $circos_cpan_dir -L $circos_cpan_dir IO::File
cpanm -l $circos_cpan_dir -L $circos_cpan_dir List::MoreUtils
cpanm -l $circos_cpan_dir -L $circos_cpan_dir List::Util
cpanm -l $circos_cpan_dir -L $circos_cpan_dir Math::Bezier
cpanm -l $circos_cpan_dir -L $circos_cpan_dir Math::BigFloat
cpanm -l $circos_cpan_dir -L $circos_cpan_dir Math::Round
cpanm -l $circos_cpan_dir -L $circos_cpan_dir Math::VecStat
cpanm -l $circos_cpan_dir -L $circos_cpan_dir Memoize
cpanm -l $circos_cpan_dir -L $circos_cpan_dir POSIX
cpanm -l $circos_cpan_dir -L $circos_cpan_dir Params::Validate
cpanm -l $circos_cpan_dir -L $circos_cpan_dir Pod::Usage
cpanm -l $circos_cpan_dir -L $circos_cpan_dir Readonly
cpanm -l $circos_cpan_dir -L $circos_cpan_dir Regexp::Common
cpanm -l $circos_cpan_dir -L $circos_cpan_dir SVG
cpanm -l $circos_cpan_dir -L $circos_cpan_dir Set::IntSpan
cpanm -l $circos_cpan_dir -L $circos_cpan_dir Statistics::Basic
cpanm -l $circos_cpan_dir -L $circos_cpan_dir Storable
cpanm -l $circos_cpan_dir -L $circos_cpan_dir Sys::Hostname
cpanm -l $circos_cpan_dir -L $circos_cpan_dir Text::Balanced
cpanm -l $circos_cpan_dir -L $circos_cpan_dir Text::Format
cpanm -l $circos_cpan_dir -L $circos_cpan_dir Time::HiRes
circos-0.69-9/bin/circos -modules

#download ExpansionHunter
cd $folder
wget https://github.com/Illumina/ExpansionHunter/releases/download/v4.0.2/ExpansionHunter-v4.0.2-linux_x86_64.tar.gz
tar xzf ExpansionHunter-v4.0.2-linux_x86_64.tar.gz
rm ExpansionHunter-v4.0.2-linux_x86_64.tar.gz

#download and build python3
mkdir -p Python3
wget https://www.python.org/ftp/python/3.6.9/Python-3.6.9.tgz
tar -zxvf Python-3.6.9.tgz
cd Python-3.6.9
./configure --prefix=$folder/Python3
make
make install
rm -R Python-3.6.9
rm Python-3.6.9.tgz

#download Splicing tools
cd $folder
spliceFolder=$folder/SplicingTools
mkdir -p $spliceFolder
cd $spliceFolder
$folder/Python3/bin/python3 -m venv splice_env
source $spliceFolder/splice_env/bin/activate
pip install cyvcf2==0.20.5 cython==0.29.21
pip install h5py==2.10.0
pip install spliceai==1.3.1
deactivate

#download REViewer
cd $folder
git clone https://github.com/Illumina/REViewer.git REViewer-0.1.1
cd REViewer-0.1.1
git fetch && git fetch --tags
git checkout v0.1.1
mkdir build
cd build
cmake ..
make
