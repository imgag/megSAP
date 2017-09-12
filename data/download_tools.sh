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
wget http://downloads.sourceforge.net/project/samtools/samtools/1.5/samtools-1.5.tar.bz2
tar xjf samtools-1.5.tar.bz2
rm samtools-1.5.tar.bz2
cd samtools-1.5
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

#download ABRA2
cd $folder
mkdir abra2-2.05
cd abra2-2.05
wget https://github.com/mozack/abra2/releases/download/v2.05/abra2-2.05.jar -O abra2.jar

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
echo "GRCh37.75.MT.codonTable : Vertebrate_Mitochondrial" >> snpEff.config
wget http://downloads.sourceforge.net/project/snpeff/databases/v4_3/snpEff_v4_3_GRCh37.75.zip
unzip snpEff_v4_3_GRCh37.75.zip
rm snpEff_v4_3_GRCh37.75.zip


