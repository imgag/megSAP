set -e
set -o pipefail
set -o verbose

folder=`pwd`/tools/
cd $folder

#download pysam, neccessary for STAR-fusion
wget https://github.com/pysam-developers/pysam/archive/v0.15.2.tar.gz
tar -xvzf v0.15.2.tar.gz
rm v0.15.2.tar.gz
cd pysam-0.15.2
export PYTHONPATH=$(pwd)/lib/python2.7/site-packages/
echo "python setup.py install --prefix=$(pwd)"
python setup.py install --prefix=$(pwd)

#download igv-reports, neccessary for STAR-fusion
wget https://github.com/igvteam/igv-reports/archive/v0.9.1.tar.gz
tar -xvzf v0.9.1.tar.gz
rm v0.9.1.tar.gz
cd igv-reports-0.9.1
export PYTHONPATH=$(pwd)/lib/python2.7/site-packages/
echo "python setup.py install --prefix=$(pwd)"
python setup.py install --prefix=$(pwd)

#download strelka2
wget https://github.com/Illumina/strelka/releases/download/v2.9.9/strelka-2.9.9.centos6_x86_64.tar.bz2
tar xjf strelka-2.9.9.centos6_x86_64.tar.bz2
rm strelka-2.9.9.centos6_x86_64.tar.bz2


#download manta
cd $folder
wget https://github.com/Illumina/manta/releases/download/v1.5.0/manta-1.5.0.centos6_x86_64.tar.bz2
tar xjf manta-1.5.0.centos6_x86_64.tar.bz2
rm manta-1.5.0.centos6_x86_64.tar.bz2
cd manta-1.5.0.centos6_x86_64
sed -i 's#referenceFasta = /illumina/development/Isis/Genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa##g' bin/configManta.py.ini


#download mantis
cd $folder
wget https://github.com/OSU-SRLab/MANTIS/archive/v1.0.4.tar.gz
tar xzf v1.0.4.tar.gz
rm v1.0.4.tar.gz
