set -e
set -o pipefail
set -o verbose

folder=`pwd`/tools/

#download RepeatMasker
cd $folder
wget http://www.repeatmasker.org/RepeatMasker-open-4-0-6.tar.gz
tar xzf RepeatMasker-open-4-0-6.tar.gz
rm -rf RepeatMasker-open-4-0-6.tar.gz

#download and build ngs-bits
cd $folder
git clone --recursive https://github.com/imgag/ngs-bits.git
cd ngs-bits
make build_3rdparty
make build_tools_release

#download and build samtools
cd $folder
wget http://downloads.sourceforge.net/project/samtools/samtools/1.3.1/samtools-1.3.1.tar.bz2
tar xjf samtools-1.3.1.tar.bz2
rm samtools-1.3.1.tar.bz2
cd samtools-1.3.1
make

#download and build freebayes
cd $folder
git clone https://github.com/ekg/freebayes.git 
cd freebayes
git checkout v1.1.0 && git submodule update --recursive --init
make

#download and build vcflib
cd $folder
git clone https://github.com/vcflib/vcflib.git
cd vcflib
git checkout v1.0.0-rc1 && git submodule update --recursive --init
make

#download ABRA
cd $folder
mkdir abra-0.97b
cd abra-0.97b
wget https://github.com/mozack/abra/releases/download/v0.97b/abra-0.97b-SNAPSHOT-jar-with-dependencies.jar -O abra.jar

#download and build samblaster
cd $folder
git clone https://github.com/GregoryFaust/samblaster.git
cd samblaster
git checkout v.0.1.24
make

#download and build BWA
cd $folder
wget http://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.15.tar.bz2
tar xjf bwa-0.7.15.tar.bz2
rm bwa-0.7.15.tar.bz2
cd bwa-0.7.15
make

#download snpEff/SnpSift
cd $folder
wget https://downloads.sourceforge.net/project/snpeff/snpEff_v4_3i_core.zip
unzip snpEff_v4_3i_core.zip
rm snpEff_v4_3i_core.zip
cd snpEff
wget http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_hg19.zip
unzip snpEff_v4_3_hg19.zip
rm snpEff_v4_3_hg19.zip


