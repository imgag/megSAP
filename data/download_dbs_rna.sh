#!/bin/bash
set -e
set -o pipefail
set -o verbose

if [ -z "$1" ]
  then
    data_folder=`pwd`
  else
    data_folder=$1
fi

src=`pwd`/../src/

#download Ensembl data
cd $data_folder/dbs/
mkdir gene_annotations
cd gene_annotations
wget -O - 'https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/__pre_July202017/GRCh37_gencode_v19_CTAT_lib_July272016.tar.gz ' | \
  gunzip -c | \
  awk '{ if($1 !~ /^#/){print "chr"$0} else{print $0} }' > GRCh37.gtf

#STAR: index genome 
cd $data_folder
mkdir -p genomes/STAR/GRCh37
$data_folder/tools/STAR_2.5.2b/bin/Linux_x86_64/STAR --runThreadN 4 --runMode genomeGenerate --genomeDir genomes/STAR/GRCh37/ --genomeFastaFiles genomes/GRCh37.fa --sjdbGTFfile dbs/gene_annotations/GRCh37.gtf

#STAR-Fusion: download and process index
cd $data_folder
mkdir -p genomes/STAR-Fusion/GRCh37
cd genomes/STAR-Fusion/GRCh37
wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh37_gencode_v19_CTAT_lib_July272016.tar.gz
tar xzf GRCh37_gencode_v19_CTAT_lib_July272016.tar.gz
rm GRCh37_gencode_v19_CTAT_lib_July272016.tar.gz
mv GRCh37_gencode_v19_CTAT_lib_July272016/* .
rmdir GRCh37_gencode_v19_CTAT_lib_July272016
export PATH=$data_folder/tools/STAR_2.5.2b/bin/Linux_x86_64/:$PATH
$data_folder/tools/STAR-Fusion-v1.0.0/FusionFilter/prep_genome_lib.pl --genome_fa ref_genome.fa --gtf ref_annot.gtf --blast_pairs blast_pairs.outfmt6.gz
                         