#!/bin/bash

root=`pwd`
folder=$root/tools/

#download and build ngs-bits
cd $folder
git clone https://github.com/imgag/ngs-bits.git
cd ngs-bits
git checkout 2020_06 && git submodule update --recursive --init
make build_3rdparty
make build_tools_release

#download and build samtools
cd $folder
wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2
tar xjf samtools-1.10.tar.bz2
rm samtools-1.10.tar.bz2
cd samtools-1.10
make

#download and build freebayes
cd $folder
git clone https://github.com/ekg/freebayes.git
cd freebayes
git checkout v1.3.2 && git submodule update --recursive --init
make

#download and build vcflib
cd $folder
git clone https://github.com/vcflib/vcflib.git
cd vcflib
git checkout v1.0.1 && git submodule update --recursive --init
make

#download ABRA2
cd $folder
mkdir abra2-2.22
cd abra2-2.22
wget https://github.com/mozack/abra2/releases/download/v2.22/abra2-2.22.jar -O abra2.jar

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

#download ClinCNV
cd $folder
git clone https://github.com/imgag/ClinCNV.git
cd ClinCNV
git fetch && git fetch --tags
git checkout 1.16.6
cd ..
mv ClinCNV ClinCNV-1.16.6

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
wget https://github.com/Illumina/ExpansionHunter/releases/download/v3.2.2/ExpansionHunter-v3.2.2-linux_x86_64.tar.gz
tar xzf ExpansionHunter-v3.2.2-linux_x86_64.tar.gz
rm ExpansionHunter-v3.2.2-linux_x86_64.tar.gz
#update variant catalog with newer github version (commit 274903d (25 Oct 2019))
wget https://raw.githubusercontent.com/Illumina/ExpansionHunter/274903d26a33cfbc546aac98c85bbfe51701fd3b/variant_catalog/grch37/variant_catalog.json -O ExpansionHunter-v3.2.2-linux_x86_64/variant_catalog/grch37/variant_catalog.json
wget https://raw.githubusercontent.com/Illumina/ExpansionHunter/274903d26a33cfbc546aac98c85bbfe51701fd3b/variant_catalog/grch38/variant_catalog.json -O ExpansionHunter-v3.2.2-linux_x86_64/variant_catalog/grch38/variant_catalog.json
wget https://raw.githubusercontent.com/Illumina/ExpansionHunter/274903d26a33cfbc546aac98c85bbfe51701fd3b/variant_catalog/hg19/variant_catalog.json -O ExpansionHunter-v3.2.2-linux_x86_64/variant_catalog/hg19/variant_catalog.json
wget https://raw.githubusercontent.com/Illumina/ExpansionHunter/274903d26a33cfbc546aac98c85bbfe51701fd3b/variant_catalog/hg38/variant_catalog.json -O ExpansionHunter-v3.2.2-linux_x86_64/variant_catalog/hg38/variant_catalog.json

