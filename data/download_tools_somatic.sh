set -e
set -o pipefail
set -o verbose

folder=`pwd`/tools/
cd $folder

#downloada strelka2
wget https://github.com/Illumina/strelka/releases/download/v2.7.1/strelka-2.7.1.centos5_x86_64.tar.bz2
tar xjf strelka-2.7.1.centos5_x86_64.tar.bz2
rm strelka-2.7.1.centos5_x86_64.tar.bz2
cd strelka-2.7.1.centos5_x86_64


#download manta
cd $folder
wget https://github.com/Illumina/manta/releases/download/v1.0.1/manta-1.0.1.centos5_x86_64.tar.bz2
tar xjf manta-1.0.1.centos5_x86_64.tar.bz2
rm manta-1.0.1.centos5_x86_64.tar.bz2
cd manta-1.0.1.centos5_x86_64
sed -i 's#referenceFasta = /illumina/development/Isis/Genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa##g' bin/configManta.py.ini