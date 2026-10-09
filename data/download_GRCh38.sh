#!/bin/bash
set -e
set -o pipefail
set -o verbose

#General description:
#  https://wiki.dnanexus.com/scientific-notes/human-genome
#Description of decoy sequences:
#  ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/README_human_reference_20110707

src=`pwd`/../src/

root=`dirname $(pwd)`

# Get data_folder path
SETTINGS_FILE=$root/settings.ini
if [ ! -f "$SETTINGS_FILE" ]; then
    SETTINGS_FILE="$root/settings.ini.default"
fi
DATA_FOLDER=$(grep -E "^data_folder" "$SETTINGS_FILE" | awk -F ' = ' '{print $2}' | sed "s|\[path\]|${root}|")

folder=$DATA_FOLDER/genomes/

# Ensure the genome folder exists
mkdir -p $folder

genome=$folder/GRCh38.fa

rm -rf $genome $genome.fai

cd $folder
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
(gunzip -c GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz |  sed -r 's/>chrM/>chrMT/g' > $genome) || true
rm GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz

php $src/Install/index_genome.php -in $genome -mask
