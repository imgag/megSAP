set -e
set -o pipefail
set -o verbose

folder=`pwd`/tools/
cd $folder

#download and build strelka
wget https://sites.google.com/site/strelkasomaticvariantcaller/home/download/strelka_workflow-1.0.15.tar.gz
tar xzf strelka_workflow-1.0.15.tar.gz
rm -rf strelka_workflow-1.0.15.tar.gz
cd strelka_workflow-1.0.15
./configure --prefix=${PWD}
make
cp etc/strelka_config_bwa_default.ini etc/strelka_config_bwa.ini
sed -i 's/isSkipDepthFilters = 0/isSkipDepthFilters = 1/g;s/maxInputDepth = 10000/maxInputDepth = -1/g' etc/strelka_config_bwa.ini
cp etc/strelka_config_bwa_default.ini etc/strelka_config_bwa_amplicon.ini
sed -i 's/ssnvNoise = 0.0000005/ssnvNoise = 0.000001/g;s/ssnvNoiseStrandBiasFrac = 0.5/ssnvNoiseStrandBiasFrac = 0/g' etc/strelka_config_bwa_amplicon.ini

#preparation for strelka2
#wget https://github.com/Illumina/strelka/releases/download/v2.7.1/strelka-2.7.1.centos5_x86_64.tar.bz2
#tar xzf strelka-2.7.1.centos5_x86_64.tar.bz2
#rm -rf strelka-2.7.1.centos5_x86_64.tar.bz2
#cd strelka-2.7.1.centos5_x86_64


#download and build manta
cd $folder
wget https://github.com/Illumina/manta/releases/download/v1.0.1/manta-1.0.1.centos5_x86_64.tar.bz2
tar xjf manta-1.0.1.centos5_x86_64.tar.bz2
rm -rf manta-1.0.1.centos5_x86_64.tar.bz2
cd manta-1.0.1.centos5_x86_64
sed -i 's#referenceFasta = /illumina/development/Isis/Genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa##g' bin/configManta.py.ini