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
git checkout 2024_02 && git submodule update --recursive --init
make build_3rdparty
make build_tools_release

#download and build python2.7 (required by manta)
cd $folder
mkdir -p Python-2.7.18
cd Python-2.7.18
wget https://www.python.org/ftp/python/2.7.18/Python-2.7.18.tgz
tar -zxvf Python-2.7.18.tgz
cd Python-2.7.18
./configure --prefix=$folder/Python-2.7.18
make
make install
cd ..
rm -R Python-2.7.18
rm Python-2.7.18.tgz
curl https://bootstrap.pypa.io/pip/2.7/get-pip.py --output get-pip.py
./bin/python2 get-pip.py
./bin/pip2 install numpy==1.16.6
./bin/pip2 install pysam==0.20.0
cd ..

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
$folder/Python-3.10.9/bin/python3 -m venv Python-3.10.9_8.2.24
source $folder/Python-3.10.9_8.2.24/bin/activate
pip install -r $root/install_deps_python.txt --require-virtualenv
deactivate
cd ..

#Python: install genome for SigProfilerExtractor - for somatic pipeline
chmod -R 777 $folder/Python-3.10.9_8.2.24/lib/python3.10/site-packages/SigProfiler*
$folder/Python-3.10.9_8.2.24/bin/python3 $root/../src/NGS/extract_signatures.py --installGenome --reference GRCh38 --in . --outFolder .


#download and build samtools
cd $folder
wget https://github.com/samtools/samtools/releases/download/1.19/samtools-1.19.tar.bz2
tar xjf samtools-1.19.tar.bz2
rm samtools-1.19.tar.bz2
cd samtools-1.19
make

#download and build bcftools
cd $folder
wget https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2
tar xjf bcftools-1.19.tar.bz2
rm bcftools-1.19.tar.bz2
cd bcftools-1.19
make

#download and build BWA
cd $folder
wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2
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

#download and build vcflib
#cd $folder
#git clone https://github.com/vcflib/vcflib.git vcflib-1.0.9
#cd vcflib-1.0.9
#git checkout v1.0.9 && git submodule update --recursive --init
#mkdir -p build && cd build
#cmake -DPYTHON_EXECUTABLE:FILEPATH=$python3_path/bin/python3 -DZIG=OFF ..
#cmake --build .
#cmake -DCMAKE_INSTALL_MANDIR=$folder/vcflib-1.0.9 --install . --prefix $folder/vcflib-1.0.9

#download and build samblaster
cd $folder
git clone https://github.com/GregoryFaust/samblaster.git samblaster-0.1.26
cd samblaster-0.1.26
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
git clone https://github.com/imgag/ClinCNV.git ClinCNV-1.18.3
cd ClinCNV-1.18.3
git checkout 1.18.3
chmod -R 777 . #if the executing user has no write permission, the error 'cannot open file Rplots.pdf' occurs
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
wget https://github.com/Illumina/interop/releases/download/v1.2.4/interop-1.2.4-Linux-GNU.tar.gz
tar xzf interop-1.2.4-Linux-GNU.tar.gz
rm interop-1.2.4-Linux-GNU.tar.gz

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

#download Splicing tools
cd $folder
spliceFolder=$folder/SplicingTools
mkdir -p $spliceFolder
cd $spliceFolder
$folder/Python-3.10.9/bin/python3 -m venv splice_env3_10
source $spliceFolder/splice_env3_10/bin/activate
pip install --upgrade pip
pip install spliceai==1.3.1
pip install tensorflow==2.11.0
deactivate
cd ..

#download REViewer
cd $folder
mkdir REViewer-v0.2.7
cd REViewer-v0.2.7
wget -O - https://github.com/Illumina/REViewer/releases/download/v0.2.7/REViewer-v0.2.7-linux_x86_64.gz | gunzip > REViewer-v0.2.7
chmod 755 REViewer-v0.2.7

#download bedtools
cd $folder
mkdir bedtools-2.31.0
cd bedtools-2.31.0
wget https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools.static
chmod 755 bedtools.static

#download illuminia ORA decompression tool
cd $folder
wget https://webdata.illumina.com/downloads/software/dragen-decompression/orad.2.6.1.tar.gz
tar xzf orad.2.6.1.tar.gz
rm orad.2.6.1.tar.gz


