#!/bin/bash
set -e
set -o pipefail
set -o verbose

root=`pwd`
tools=$root/tools/
dbs=$root/dbs/

vep_install_dir=$tools/ensembl-vep-release-93.2/
vep_cpan_dir=$vep_install_dir/cpan/
vep_data_dir=$dbs/ensembl-vep-93/

# download ensembl-vep
cd $tools
wget https://github.com/Ensembl/ensembl-vep/archive/release/93.2.tar.gz
tar xzf 93.2.tar.gz

# install PERL dependencies
mkdir -p $vep_cpan_dir
cpanm -l $vep_cpan_dir -L $vep_cpan_dir Set::IntervalTree PerlIO::gzip DBI

# install BigFile support (for BigWig support needed to annotate phyloP)
cd $vep_install_dir
export KENT_SRC=$PWD/kent-335_base/src
export MACHTYPE=$(uname -m)
export CFLAGS="-fPIC"
wget https://github.com/ucscGenomeBrowser/kent/archive/v335_base.tar.gz
tar xzf v335_base.tar.gz
cd $KENT_SRC/lib
echo 'CFLAGS="-fPIC"' > ../inc/localEnvironment.mk
make clean && make
cd ../jkOwnLib
make clean && make
cpanm -l $vep_cpan_dir -L $vep_cpan_dir Bio::DB::BigFile

#download reference genome for VEP //@TODO can we perhaps use the genome from the genomes folder? Is re-zipping with bgzip really necessary?
cd $vep_data_dir
mkdir -p fasta
cd fasta
curl -O ftp://ftp.ensembl.org/pub/grch37/release-93/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
gzip -d Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
bgzip Homo_sapiens.GRCh37.dna.primary_assembly.fa

#download VEP cache data
cd $vep_data_dir
mkdir -p ftp
cd ftp
wget ftp://ftp.ensembl.org/pub/release-93/variation/VEP/homo_sapiens_vep_93_GRCh37.tar.gz

#install ensembl-vep
cd $vep_install_dir
perl INSTALL.pl --SPECIES homo_sapiens --ASSEMBLY GRCh37 --AUTO acp --PLUGINS REVEL,FATHMM_MKL,CADD,dbscSNV,GeneSplicer,MaxEntScan --NO_UPDATE --CACHEDIR $vep_data_dir/cache --CACHEURL $vep_data_dir/ftp --FASTAURL $vep_data_dir/fasta
cp $vep_data_dir/cache/Plugins/*.pm $vep_install_dir/modules/ #should not be necessary - probably a bug in the VEP installation script when using the CACHEDIR option (MS)

# install MaxEntScan (for MaxEntScan plugin)
cd $vep_install_dir
mkdir -p MaxEntScan
cd MaxEntScan
wget http://genes.mit.edu/burgelab/maxent/download/fordownload.tar.gz
tar xzf fordownload.tar.gz
mv fordownload/* .
rm -rf fordownload*

# install GeneSplicer (for GeneSplicer plugin)
cd $vep_install_dir
wget ftp://ccb.jhu.edu/pub/software/genesplicer/GeneSplicer.tar.gz
tar -xzf GeneSplicer.tar.gz
rm -rf GeneSplicer.tar.gz
