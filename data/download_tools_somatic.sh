set -e
set -o pipefail
set -o verbose

folder=`pwd`/tools/
cd $folder

#downloada strelka2
wget https://github.com/Illumina/strelka/releases/download/v2.9.9/strelka-2.9.9.centos6_x86_64.tar.bz2
tar xjf strelka-2.9.9.centos6_x86_64.tar.bz2
rm strelka-2.9.9.centos6_x86_64.tar.bz2


#download manta
cd $folder
wget https://github.com/Illumina/manta/releases/download/v1.4.0/manta-1.4.0.centos6_x86_64.tar.bz2
tar xjf manta-1.4.0.centos6_x86_64.tar.bz2
rm manta-1.4.0.centos6_x86_64.tar.bz2
cd manta-1.4.0.centos6_x86_64
sed -i 's#referenceFasta = /illumina/development/Isis/Genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa##g' bin/configManta.py.ini


#download mantis
cd $folder
wget https://github.com/OSU-SRLab/MANTIS/archive/v1.0.4.tar.gz
tar xzf v1.0.4.tar.gz
rm v1.0.4.tar.gz