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

#download Ensembl data
cd $data_folder
mkdir -p dbs/gene_annotations
wget -O - 'ftp://ftp.ensembl.org/pub/grch37/update/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gtf.gz' | \
  gzip -cd | \
  awk '{ if ($1 !~ /^#/) { print "chr"$0 } else { print $0 } }' > dbs/gene_annotations/GRCh37.gtf

#STAR: index genome
cd $data_folder
mkdir -p genomes/STAR/GRCh37
$data_folder/tools/STAR-2.6.1a/bin/Linux_x86_64/STAR --runThreadN 4 --runMode genomeGenerate --genomeDir genomes/STAR/GRCh37/ --genomeFastaFiles genomes/GRCh37.fa --sjdbGTFfile dbs/gene_annotations/GRCh37.gtf

#STAR-Fusion: download pre-computed index
cd $data_folder
mkdir -p genomes/STAR-Fusion
cd genomes/STAR-Fusion
wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh37_v19_CTAT_lib_Feb092018.source_data.tar.gz
tar xzf GRCh37_v19_CTAT_lib_Feb092018.source_data.tar.gz
cd GRCh37_v19_CTAT_lib_Feb092018
export PATH=$data_folder/tools/STAR-2.6.1c/bin/Linux_x86_64/:$PATH
$data_folder/tools/STAR-Fusion-v1.5.0/FusionFilter/prep_genome_lib.pl \
  --genome_fa ref_genome.fa \
  --gtf ref_annot.gtf \
  --fusion_annot_lib CTAT_HumanFusionLib.v0.1.0.dat.gz \
  --annot_filter_rule AnnotFilterRule.pm \
  --pfam_db PFAM.domtblout.dat.gz
cd ..
mv GRCh37_v19_CTAT_lib_Feb092018/ctat_genome_lib_build_dir GRCh37
rm -r GRCh37_v19_CTAT_lib_Feb092018 GRCh37_v19_CTAT_lib_Feb092018.source_data.tar.gz
