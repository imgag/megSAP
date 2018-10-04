#!/bin/bash

root=`pwd`

#download reference genome for VEP //@TODO can we perhaps use the genome from the genomes folder?
cd $root
rm -rf fasta
mkdir fasta
cd fasta
curl -O ftp://ftp.ensembl.org/pub/grch37/release-93/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
gzip -d Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
bgzip Homo_sapiens.GRCh37.dna.primary_assembly.fa

#download VEP cache data
cd $root
rm -rf ftp
mkdir ftp
cd ftp
wget ftp://ftp.ensembl.org/pub/release-93/variation/VEP/homo_sapiens_vep_93_GRCh37.tar.gz
tar xzf homo_sapiens_vep_93_GRCh37.tar.gz
