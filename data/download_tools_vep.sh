#!/bin/bash
set -e
set -o pipefail
set -o verbose

mkdir -p `pwd`/tools
folder=`pwd`/tools/

cd $folder

# download ensembl-vep
git clone https://github.com/Ensembl/ensembl-vep.git --depth 1
cd ensembl-vep

# install ensembl dependencies
cpanm DBI DBD::mysql

# build ensembl-vep without cache or plugins
perl INSTALL.pl -a a

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

# install BigWig
mkdir -p $HOME/cpanm
export PERL5LIB=$PERL5LIB:$HOME/cpanm/lib/perl5
cpanm -l $HOME/cpanm Bio::DB::BigFile
