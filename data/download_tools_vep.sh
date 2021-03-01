#!/bin/bash
set -e
set -o pipefail
set -o verbose

root=`pwd`
tools=$root/tools/
dbs=$root/dbs/

#Ignore this - used for local installation
tools=/mnt/share/opt/
# tools=/home/bioinf/
dbs=/mnt/share/data/dbs/

vep_install_dir=$tools/ensembl-vep-release-103.1/
vep_cpan_dir=$vep_install_dir/cpan/
vep_data_dir=$dbs/ensembl-vep-103/

# download ensembl-vep
cd $tools
wget https://github.com/Ensembl/ensembl-vep/archive/release/103.1.tar.gz
mkdir $vep_install_dir
tar -C $vep_install_dir --strip-components=1 -xzf 103.1.tar.gz
rm 103.1.tar.gz

#install dependencies
mkdir -p $vep_cpan_dir
cpanm -l $vep_cpan_dir -L $vep_cpan_dir Set::IntervalTree URI::Escape DB_File Carp::Assert JSON::XS PerlIO::gzip DBI

#install BigWig support (needed to annotate phyloP)
cd $vep_install_dir
export KENT_SRC=$vep_install_dir/kent-335_base/src
export MACHTYPE=$(uname -m)
export CFLAGS="-fPIC"
wget https://github.com/ucscGenomeBrowser/kent/archive/v335_base.tar.gz
tar xzf v335_base.tar.gz
rm v335_base.tar.gz
cd $KENT_SRC/lib
echo 'CFLAGS="-fPIC"' > $KENT_SRC/inc/localEnvironment.mk
make clean && make
cd $KENT_SRC/jkOwnLib
make clean && make
cpanm -l $vep_cpan_dir -L $vep_cpan_dir Bio::DB::BigFile

#download VEP cache data
mkdir -p $vep_data_dir
cd $vep_data_dir
mkdir -p ftp
cd ftp
wget ftp://ftp.ensembl.org/pub/release-103/variation/indexed_vep_cache/homo_sapiens_vep_103_GRCh37.tar.gz
wget ftp://ftp.ensembl.org/pub/release-103/variation/indexed_vep_cache/homo_sapiens_refseq_vep_103_GRCh37.tar.gz

#install ensembl-vep
PERL5LIB=$vep_install_dir/Bio/:$vep_cpan_dir/lib/perl5/:$PERL5LIB
cd $vep_install_dir
perl INSTALL.pl --SPECIES homo_sapiens,homo_sapiens_refseq --ASSEMBLY GRCh37 --AUTO acp --PLUGINS REVEL,FATHMM_MKL,CADD,dbscSNV,MaxEntScan --NO_UPDATE --NO_BIOPERL --CACHEDIR $vep_data_dir/cache --CACHEURL $vep_data_dir/ftp --NO_TEST
cp $vep_data_dir/cache/Plugins/*.pm $vep_install_dir/modules/ #should not be necessary - probably a bug in the VEP installation script when using the CACHEDIR option (MS)

# install MaxEntScan (for MaxEntScan plugin)
cd $vep_install_dir
mkdir -p MaxEntScan
cd MaxEntScan
wget http://hollywood.mit.edu/burgelab/maxent/download/fordownload.tar.gz
tar xzf fordownload.tar.gz
mv fordownload/* .
rm -rf fordownload*
chmod -R 755 $vep_install_dir/MaxEntScan
