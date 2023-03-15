#!/bin/bash
set -e
set -o pipefail
set -o verbose

root=`pwd`
folder=$root/tools/
dbs=$root/dbs/

#Ignore this - used for local installation
#folder=/mnt/storage2/megSAP/tools/
#dbs=/mnt/storage2/megSAP/data/dbs

vep_install_dir=$folder/ensembl-vep-release-107.0/
vep_data_dir=$dbs/ensembl-vep-107/
cpan_dir=$folder/perl_cpan/

# download ensembl-vep
cd $folder
wget https://github.com/Ensembl/ensembl-vep/archive/release/107.0.tar.gz
mkdir $vep_install_dir
tar -C $vep_install_dir --strip-components=1 -xzf 107.0.tar.gz
rm 107.0.tar.gz

#install dependencies
mkdir -p $cpan_dir
cpanm -l $cpan_dir -L $cpan_dir Set::IntervalTree URI::Escape DB_File Carp::Assert JSON::XS PerlIO::gzip DBI

#download VEP cache data
mkdir -p $vep_data_dir
cd $vep_data_dir
mkdir -p ftp
cd ftp
wget ftp://ftp.ensembl.org/pub/release-107/variation/indexed_vep_cache/homo_sapiens_vep_107_GRCh38.tar.gz

#install ensembl-vep
PERL5LIB=$vep_install_dir/Bio/:$cpan_dir/lib/perl5/:$PERL5LIB
cd $vep_install_dir
perl INSTALL.pl --SPECIES homo_sapiens --ASSEMBLY GRCh38 --AUTO acp --PLUGINS REVEL,CADD,MaxEntScan --NO_UPDATE --NO_BIOPERL --CACHEDIR $vep_data_dir/cache --CACHEURL $vep_data_dir/ftp --NO_TEST --NO_HTSLIB
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
