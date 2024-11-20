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
wget -O - https://www.proteinatlas.org/download/rna_tissue_consensus.tsv.zip | gunzip > rna_tissue_consensus_v22.tsv

#download Ensembl data in GTF format - KEEP AT ENSEMBL VERSION 109, DB import for RNA works on Transcript base and will break if the transcripts change.
cd $data_folder
mkdir -p dbs/gene_annotations
cd dbs/gene_annotations
wget -O - 'https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz' | gzip -cd | awk '{ if ($$1 !~ /^#/) { print "chr"$0 } else { print $0 } }' > GRCh38.gtf

#create hemoglobin FASTA file
cd $data_folder
cd misc
wget -O - 'https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz' | zcat | \
	awk -v RS=">" -v FS="\n" '$$1 ~ / gene_symbol:(HBA1|HBA2|HBB) / { print ">"$$1; {for (i=2; i<=NF; i++) printf("%s", $$i)}; printf("\n") }' | \
	sed '/^>/s/ /|kraken:taxid|9606 /' \
	> human_hemoglobin_tx.fa
