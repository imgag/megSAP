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
	genome=`pwd`/genomes/GRCh37.fa
  else
	mkdir -p $1/genomes/
	genome=$1/genomes/GRCh37.fa
fi
rm -rf $genome $genome.fai

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
(gunzip -c hs37d5.fa.gz > $genome) || true
rm hs37d5.fa.gz

php $src/Tools/index_genome.php -in $genome
