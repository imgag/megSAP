set -e
set -o pipefail
set -o verbose

folder=`pwd`/tools/

#download STAR
cd $folder
git clone https://github.com/alexdobin/STAR.git STAR-2.5.4a
cd STAR-2.5.4a
git checkout 2.5.4a
git submodule update --recursive --init

#download STAR-Fusion
cd $folder
wget https://github.com/STAR-Fusion/STAR-Fusion/releases/download/STAR-Fusion-v1.2.0/STAR-Fusion-v1.2.0.FULL.tar.gz
tar xzf STAR-Fusion-v1.2.0.FULL.tar.gz
rm STAR-Fusion-v1.2.0.FULL.tar.gz

#download subread (featureCounts, subjunc)
cd $folder
wget http://downloads.sourceforge.net/project/subread/subread-1.6.0/subread-1.6.0-Linux-x86_64.tar.gz
tar xzf subread-1.6.0-Linux-x86_64.tar.gz
rm subread-1.6.0-Linux-x86_64.tar.gz
mv subread-1.6.0-Linux-x86_64 subread-1.6.0

#download skewer (single-end adapter trimming)
cd $folder
git clone https://github.com/relipmoc/skewer.git skewer_0.2.2
cd skewer_0.2.2
git checkout 0.2.2
make

#download salmon
wget https://github.com/COMBINE-lab/salmon/releases/download/v0.9.1/Salmon-0.9.1_linux_x86_64.tar.gz
tar xzf Salmon-0.9.1_linux_x86_64.tar.gz
rm Salmon-0.9.1_linux_x86_64.tar.gz
mv Salmon-latest_linux_x86_64 salmon-0.9.1