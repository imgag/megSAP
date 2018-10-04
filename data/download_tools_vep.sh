#!/bin/bash
folder=`pwd`/tools/

if [ `basename "$PWD"` != "tools" ]; then
  cd $folder
fi

# download ensembl-vep
git clone https://github.com/Ensembl/ensembl-vep.git --depth 1
cd ensembl-vep
git checkout release/93

# install ensembl dependencies
cpanm DBI DBD::mysql

# build ensembl-vep without cache or plugins //@TODO here we need to provice the ftp/cache parameter
perl INSTALL.pl -a ap -g dbscSNV,GeneSplicer,MaxEntScan,REVEL,FATHMM_MKL,CADD

# build BigWig
cd $folder
export KENT_SRC=$PWD/kent-335_base/src
export MACHTYPE=$(uname -m)
export CFLAGS="-fPIC"
export MYSQLINC=`mysql_config --include | sed -e 's/^-I//g'`
export MYSQLLIBS=`mysql_config --libs`

echo 'Getting jksrc'
if [ ! -f v335_base.tar.gz ]; then
  wget https://github.com/ucscGenomeBrowser/kent/archive/v335_base.tar.gz
  tar xzf v335_base.tar.gz
fi

cd $KENT_SRC/lib
echo 'CFLAGS="-fPIC"' > ../inc/localEnvironment.mk

make clean && make
cd ../jkOwnLib
make clean && make

# install BigWig support (for phyloP custom annotation)
mkdir -p $HOME/cpanm
export PERL5LIB=$PERL5LIB:$HOME/cpanm/lib/perl5
cpanm -l $HOME/cpanm Bio::DB::BigFile

# install MaxEntScan (for MaxEntScan plugin)
cd $folder/ensembl-vep
mkdir -p MaxEntScan
cd MaxEntScan
wget http://genes.mit.edu/burgelab/maxent/download/fordownload.tar.gz
tar xzf fordownload.tar.gz
mv fordownload/* .
rm -rf fordownload*

# install GeneSplicer (for GeneSplicer plugin)
cd $folder/ensembl-vep
wget ftp://ccb.jhu.edu/pub/software/genesplicer/GeneSplicer.tar.gz
tar -xzf GeneSplicer.tar.gz
rm -rf GeneSplicer.tar.gz