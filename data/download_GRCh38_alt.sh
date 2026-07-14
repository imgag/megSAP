#!/bin/bash
set -e
set -o pipefail
set -o verbose

#General description:
#  https://wiki.dnanexus.com/scientific-notes/human-genome
#Description of decoy sequences:
#  ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/README_human_reference_20110707

src=`pwd`/../src/

if [ -z "$1" ]
  then
	genome=`pwd`/genomes/GRCh38_alt.fa
  else
	mkdir -p $1/genomes/
	genome=$1/genomes/GRCh38_alt.fa
fi
rm -rf $genome $genome.fai

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz
(gunzip -c GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz |  sed -r 's/>chrM/>chrMT/g' > $genome) || true
rm GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz

php $src/Install/index_genome.php -in $genome -mask
