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
wget https://www.proteinatlas.org/download/rna_tissue_hpa.tsv.zip
unzip -o rna_tissue_hpa.tsv.zip
rm rna_tissue_hpa.tsv.zip

#download Ensembl data
cd $data_folder
mkdir -p dbs/gene_annotations
wget -O - 'http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz' | \
  gzip -cd | \
  awk '{ if ($1 !~ /^#/) { print "chr"$0 } else { print $0 } }' > dbs/gene_annotations/GRCh38.gtf

#STAR: index genome
cd $data_folder
mkdir -p genomes/STAR/GRCh38
$data_folder/tools/STAR-2.7.9a/bin/Linux_x86_64/STAR \
--runThreadN 20 \
--runMode genomeGenerate \
--genomeDir genomes/STAR/GRCh38/ \
--genomeFastaFiles genomes/GRCh38.fa \
--sjdbGTFfile dbs/gene_annotations/GRCh38.gtf
