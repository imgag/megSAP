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

#download reference data for gene expression
cd $data_folder
mkdir -p dbs/gene_expression
cd dbs/gene_expression
#change version number on update
wget -O - https://www.proteinatlas.org/download/rna_tissue_hpa.tsv.zip | gunzip > rna_tissue_hpa_v21.1.tsv

#download Ensembl data
cd $data_folder
mkdir -p dbs/Ensembl
wget -O - 'http://ftp.ensembl.org/pub/release-107/gtf/homo_sapiens/Homo_sapiens.GRCh38.107.gtf.gz' | \
  gzip -cd | \
  awk '{ if ($1 !~ /^#/) { print "chr"$0 } else { print $0 } }' > dbs/Ensembl/Homo_sapiens.GRCh38.107.chr.gtf

#STAR: index genome
cd $data_folder
mkdir -p genomes/STAR/GRCh38
$data_folder/tools/STAR-2.7.10a/bin/Linux_x86_64/STAR \
--runThreadN 20 \
--runMode genomeGenerate \
--genomeDir genomes/STAR/GRCh38/ \
--genomeFastaFiles genomes/GRCh38.fa \
--sjdbGTFfile dbs/Ensembl/Homo_sapiens.GRCh38.107.chr.gtf
