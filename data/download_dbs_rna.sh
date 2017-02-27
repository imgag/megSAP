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

#download UCSC data
cd $data_folder
mkdir UCSC
cd UCSC
wget -O - http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz | gunzip > refGene.txt
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToGtf
chmod 775 genePredToGtf
cut -f2- refGene.txt | ./genePredToGtf file stdin refGene.gtf
rm genePredToGtf

#STAR: index genome 
cd $data_folder
mkdir -p genomes/STAR/hg19
tools/STAR_2.5.2b/bin/Linux_x86_64/STAR --runThreadN 4 --runMode genomeGenerate --genomeDir genomes/STAR/hg19/ --genomeFastaFiles genomes/hg19.fa --sjdbGTFfile dbs/UCSC/refGene.gtf

#STAR-Fusion: download and process index
cd $data_folder
mkdir -p genomes/STAR-Fusion/hg19
cd genomes/STAR-Fusion/hg19
wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh37_gencode_v19_CTAT_lib_July272016.tar.gz
tar xzf GRCh37_gencode_v19_CTAT_lib_July272016.tar.gz
rm GRCh37_gencode_v19_CTAT_lib_July272016.tar.gz
mv GRCh37_gencode_v19_CTAT_lib_July272016/* .
rmdir GRCh37_gencode_v19_CTAT_lib_July272016
export PATH=$data_folder/tools/STAR_2.5.2b/bin/Linux_x86_64/:$PATH
$data_folder/tools/STAR-Fusion-v1.0.0/FusionFilter/prep_genome_lib.pl --genome_fa ref_genome.fa --gtf ref_annot.gtf --blast_pairs blast_pairs.outfmt6.gz
                         