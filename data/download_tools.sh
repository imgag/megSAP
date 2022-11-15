#!/bin/bash

root=`pwd`
folder=$root/tools/
cpan_dir=$folder/perl_cpan/

#download and build ngs-bits
cd $folder
git clone https://github.com/imgag/ngs-bits.git
cd ngs-bits
git checkout 2022_10 && git submodule update --recursive --init
make build_3rdparty
make build_tools_release

#download and build samtools
cd $folder
wget https://github.com/samtools/samtools/releases/download/1.15.1/samtools-1.15.1.tar.bz2
tar xjf samtools-1.15.1.tar.bz2
rm samtools-1.15.1.tar.bz2
cd samtools-1.15.1
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
wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2
tar xjf bwa-mem2-2.2.1_x64-linux.tar.bz2
rm bwa-mem2-2.2.1_x64-linux.tar.bz2


#download ABRA2
cd $folder
mkdir abra2-2.23
cd abra2-2.23
wget https://github.com/mozack/abra2/releases/download/v2.23/abra2-2.23.jar -O abra2.jar

#download freebayes
cd $folder
mkdir freebayes-1.3.6
cd freebayes-1.3.6
wget -O - https://github.com/freebayes/freebayes/releases/download/v1.3.6/freebayes-1.3.6-linux-amd64-static.gz | gunzip > freebayes
chmod 755 freebayes

#download and build vcflib
cd $folder
git clone https://github.com/vcflib/vcflib.git vcflib-1.0.3
cd vcflib-1.0.3
git checkout v1.0.3 && git submodule update --recursive --init
mkdir -p build && cd build
cmake ..
cmake --build .
cmake --install .

#download and build samblaster
cd $folder
git clone https://github.com/GregoryFaust/samblaster.git
cd samblaster
git checkout v.0.1.26
make

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
git clone https://github.com/imgag/ClinCNV.git ClinCNV-1.18.0
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
mkdir -p $cpan_dir
cpanm -l $cpan_dir -L $cpan_dir Carp Clone Config::General Cwd Data::Dumper Digest::MD5 File::Basename File::Spec::Functions File::Temp FindBin Font::TTF::Font GD GD::Polyline Getopt::Long IO::File List::MoreUtils List::Util Math::Bezier Math::BigFloat Math::Round Math::VecStat Memoize POSIX Params::Validate Pod::Usage Readonly Regexp::Common SVG Set::IntSpan Statistics::Basic Storable Sys::Hostname Text::Balanced Text::Format Time::HiRes
circos-0.69-9/bin/circos -modules

#download ExpansionHunter
cd $folder
wget https://github.com/Illumina/ExpansionHunter/releases/download/v5.0.0/ExpansionHunter-v5.0.0-linux_x86_64.tar.gz
tar xzf ExpansionHunter-v5.0.0-linux_x86_64.tar.gz
rm ExpansionHunter-v5.0.0-linux_x86_64.tar.gz

#download and build python3
cd $folder
mkdir -p Python3
wget https://www.python.org/ftp/python/3.6.9/Python-3.6.9.tgz
tar -zxvf Python-3.6.9.tgz
cd Python-3.6.9
./configure --prefix=$folder/Python3
make
make install
cd ..
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
mkdir REViewer-v0.2.7
cd REViewer-v0.2.7
wget -O - https://github.com/Illumina/REViewer/releases/download/v0.2.7/REViewer-v0.2.7-linux_x86_64.gz | gunzip > REViewer-v0.2.7
chmod 755 REViewer-v0.2.7

#download bedtools
cd $folder
mkdir bedtools-2.30.0
cd bedtools-2.30.0
wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary
chmod 755 bedtools.static.binary
