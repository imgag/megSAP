set -e
set -o pipefail
set -o verbose

folder=`pwd`/tools/

#TODO tools for indel realignment?: GATK
#TODO tools for duplicates?: picard

#download STAR
cd $folder
git clone https://github.com/alexdobin/STAR.git
mv STAR STAR_2.5.2b
cd STAR_2.5.2b
git checkout 2.5.2b && git submodule update --recursive --init

#download STAR-Fusion
cd $folder
wget https://github.com/STAR-Fusion/STAR-Fusion/releases/download/v1.0.0/STAR-Fusion-v1.0.0.FULL.tar.gz
tar xzf STAR-Fusion-v1.0.0.FULL.tar.gz
rm -rf STAR-Fusion-v1.0.0.FULL.tar.gz


#download featureCounts
cd $folder
wget http://downloads.sourceforge.net/project/subread/subread-1.5.1/subread-1.5.1-Linux-x86_64.tar.gz
tar xzf subread-1.5.1-Linux-x86_64.tar.gz
rm -rf subread-1.5.1-Linux-x86_64.tar.gz
mv subread-1.5.1-Linux-x86_64 subread-1.5.1

#download skewer (single-end adapter trimming)
cd $folder
git clone https://github.com/relipmoc/skewer.git
mv skewer skewer_0.2.2
cd skewer_0.2.2
git checkout 0.2.2
make

