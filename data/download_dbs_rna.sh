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
wget -O - 'ftp://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.chr_patch_hapl_scaff.gtf.gz' | \
  gzip -cd | \
  awk '{ if ($1 !~ /^#/) { print "chr"$0 } else { print $0 } }' > dbs/gene_annotations/GRCh37.gtf

#STAR: index genome
cd $data_folder
mkdir -p genomes/STAR/GRCh37
$data_folder/tools/STAR-2.7.9a/bin/Linux_x86_64/STAR \
--runThreadN 20 \
--runMode genomeGenerate \
--genomeDir genomes/STAR/GRCh37/ \
--genomeFastaFiles genomes/GRCh37.fa \
--sjdbGTFfile dbs/gene_annotations/GRCh37.gtf

#STAR-Fusion: download pre-computed index
cd $data_folder
mkdir -p genomes/STAR-Fusion
cd genomes/STAR-Fusion
wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/__genome_libs_StarFv1.9/GRCh37_gencode_v19_CTAT_lib_Apr032020.plug-n-play.tar.gz
tar xzf GRCh37_gencode_v19_CTAT_lib_Apr032020.plug-n-play.tar.gz
mv GRCh37_gencode_v19_CTAT_lib_Apr032020.plug-n-play/ctat_genome_lib_build_dir GRCh37
rm -r GRCh37_gencode_v19_CTAT_lib_Apr032020.plug-n-play GRCh37_gencode_v19_CTAT_lib_Apr032020.plug-n-play.tar.gz
