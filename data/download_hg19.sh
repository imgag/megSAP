#!/bin/bash
set -e
set -o pipefail
set -o verbose

if [ -z "$1" ]
  then
	genome=`pwd`/genomes/hg19.fa
  else
	mkdir -p $1/genomes/
	genome=$1/genomes/hg19.fa
fi
rm -rf $genome $genome.fai

path="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes"

wget -O - $path/chrM.fa.gz | gunzip  >> $genome
wget -O - $path/chr1.fa.gz | gunzip  >> $genome
wget -O - $path/chr2.fa.gz | gunzip  >> $genome
wget -O - $path/chr3.fa.gz | gunzip  >> $genome
wget -O - $path/chr4.fa.gz | gunzip  >> $genome
wget -O - $path/chr5.fa.gz | gunzip  >> $genome
wget -O - $path/chr6.fa.gz | gunzip  >> $genome
wget -O - $path/chr7.fa.gz | gunzip  >> $genome
wget -O - $path/chr8.fa.gz | gunzip  >> $genome
wget -O - $path/chr9.fa.gz | gunzip  >> $genome
wget -O - $path/chr10.fa.gz | gunzip  >> $genome
wget -O - $path/chr11.fa.gz | gunzip  >> $genome
wget -O - $path/chr12.fa.gz | gunzip  >> $genome
wget -O - $path/chr13.fa.gz | gunzip  >> $genome
wget -O - $path/chr14.fa.gz | gunzip  >> $genome
wget -O - $path/chr15.fa.gz | gunzip  >> $genome
wget -O - $path/chr16.fa.gz | gunzip  >> $genome
wget -O - $path/chr17.fa.gz | gunzip  >> $genome
wget -O - $path/chr18.fa.gz | gunzip  >> $genome
wget -O - $path/chr19.fa.gz | gunzip  >> $genome
wget -O - $path/chr20.fa.gz | gunzip  >> $genome
wget -O - $path/chr21.fa.gz | gunzip  >> $genome
wget -O - $path/chr22.fa.gz | gunzip  >> $genome
wget -O - $path/chrX.fa.gz | gunzip  >> $genome
wget -O - $path/chrY.fa.gz | gunzip  >> $genome
wget -O - $path/chr1_gl000191_random.fa.gz | gunzip  >> $genome
wget -O - $path/chr1_gl000192_random.fa.gz | gunzip  >> $genome
#wget -O - $path/chr4_ctg9_hap1.fa.gz | gunzip  >> $genome
wget -O - $path/chr4_gl000193_random.fa.gz | gunzip  >> $genome
wget -O - $path/chr4_gl000194_random.fa.gz | gunzip  >> $genome
#wget -O - $path/chr6_apd_hap1.fa.gz | gunzip  >> $genome
#wget -O - $path/chr6_cox_hap2.fa.gz | gunzip  >> $genome
#wget -O - $path/chr6_dbb_hap3.fa.gz | gunzip  >> $genome
#wget -O - $path/chr6_mann_hap4.fa.gz | gunzip  >> $genome
#wget -O - $path/chr6_mcf_hap5.fa.gz | gunzip  >> $genome
#wget -O - $path/chr6_qbl_hap6.fa.gz | gunzip  >> $genome
#wget -O - $path/chr6_ssto_hap7.fa.gz | gunzip  >> $genome
wget -O - $path/chr7_gl000195_random.fa.gz | gunzip  >> $genome
wget -O - $path/chr8_gl000196_random.fa.gz | gunzip  >> $genome
wget -O - $path/chr8_gl000197_random.fa.gz | gunzip  >> $genome
wget -O - $path/chr9_gl000198_random.fa.gz | gunzip  >> $genome
wget -O - $path/chr9_gl000199_random.fa.gz | gunzip  >> $genome
wget -O - $path/chr9_gl000200_random.fa.gz | gunzip  >> $genome
wget -O - $path/chr9_gl000201_random.fa.gz | gunzip  >> $genome
wget -O - $path/chr11_gl000202_random.fa.gz | gunzip  >> $genome
#wget -O - $path/chr17_ctg5_hap1.fa.gz | gunzip  >> $genome
wget -O - $path/chr17_gl000203_random.fa.gz | gunzip  >> $genome
wget -O - $path/chr17_gl000204_random.fa.gz | gunzip  >> $genome
wget -O - $path/chr17_gl000205_random.fa.gz | gunzip  >> $genome
wget -O - $path/chr17_gl000206_random.fa.gz | gunzip  >> $genome
wget -O - $path/chr18_gl000207_random.fa.gz | gunzip  >> $genome
wget -O - $path/chr19_gl000208_random.fa.gz | gunzip  >> $genome
wget -O - $path/chr19_gl000209_random.fa.gz | gunzip  >> $genome
wget -O - $path/chr21_gl000210_random.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000211.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000212.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000213.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000214.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000215.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000216.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000217.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000218.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000219.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000220.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000221.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000222.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000223.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000224.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000225.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000226.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000227.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000228.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000229.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000230.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000231.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000232.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000233.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000234.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000235.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000236.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000237.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000238.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000239.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000240.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000241.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000242.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000243.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000244.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000245.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000246.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000247.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000248.fa.gz | gunzip  >> $genome
wget -O - $path/chrUn_gl000249.fa.gz | gunzip  >> $genome

php `pwd`/../src/Tools/index_genome.php -in $genome
