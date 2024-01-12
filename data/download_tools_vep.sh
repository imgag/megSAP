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

vep_install_dir=$folder/ensembl-vep-release-110.1/
vep_data_dir=$dbs/ensembl-vep-110/
cpan_dir=$folder/perl_cpan/

# download ensembl-vep
cd $folder
wget https://github.com/Ensembl/ensembl-vep/archive/release/110.1.tar.gz
mkdir $vep_install_dir
tar -C $vep_install_dir --strip-components=1 -xzf 110.1.tar.gz
rm 110.1.tar.gz

#install dependencies
mkdir -p $cpan_dir
cpanm -l $cpan_dir -L $cpan_dir Set::IntervalTree URI::Escape DB_File Carp::Assert JSON::XS PerlIO::gzip DBI

#download VEP cache data
mkdir -p $vep_data_dir
cd $vep_data_dir
mkdir -p ftp
cd ftp
wget https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_vep_110_GRCh38.tar.gz

#install ensembl-vep
PERL5LIB=$vep_install_dir/Bio/:$cpan_dir/lib/perl5/:$PERL5LIB
cd $vep_install_dir
perl INSTALL.pl --SPECIES homo_sapiens --ASSEMBLY GRCh38 --AUTO ac --NO_UPDATE --NO_BIOPERL --CACHEDIR $vep_data_dir/cache --CACHEURL $vep_data_dir/ftp --NO_TEST --NO_HTSLIB
