#!/bin/bash
set -e
set -o pipefail
set -o verbose

if [ -z "$1" ]
  then
	dbs=`pwd`/dbs/
  else	
	dbs=$1/dbs/
	mkdir -p $dbs
fi


src=`pwd`/../src/
tools=`pwd`/tools/
ngsbits=$tools/ngs-bits/bin
vcflib=$tools/vcflib/bin
genome=`pwd`/genomes/GRCh37.fa

#Install REPEATMASKER - http://www.repeatmasker.org/species/hg.html
cd $dbs
mkdir RepeatMasker
cd RepeatMasker
wget -O - http://www.repeatmasker.org/genomes/hg19/RepeatMasker-rm405-db20140131/hg19.fa.out.gz | gunzip > hg19.fa.out
perl $tools/RepeatMasker/util/rmOutToGFF3.pl hg19.fa.out > RepeatMasker.gff
cat RepeatMasker.gff | php $src/Tools/db_converter_repeatmasker.php | $ngsbits/BedSort | bgzip > RepeatMasker.bed
tabix -p bed RepeatMasker.bed.gz
rm -rf hg19.fa.out RepeatMasker.gff
cd_st
#Install ClinVar - https://www.ncbi.nlm.nih.gov/clinvar/
cd $dbs
mkdir ClinVar
cd ClinVar
wget -O - ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/archive_2.0/2018/clinvar_20180701.vcf.gz | gunzip | php $src/Tools/db_converter_clinvar.php | bgzip > clinvar_20180701_converted.vcf.gz
tabix -p vcf clinvar_20180701_converted.vcf.gz

#Install gnomAD genome data - https://www.ensembl.org/info/docs/tools/vep/script/vep_example.html#gnomad
cd $dbs
mkdir gnomAD
cd gnomAD
wget ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh37/variation_genotype/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz
wget ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh37/variation_genotype/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz.tbi

#Install phyloP - https://www.ensembl.org/info/docs/tools/vep/script/vep_example.html#gerp
cd $dbs
mkdir phyloP
cd phyloP
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/phyloP100way/hg19.100way.phyloP100way.bw

#Install CADD - http://cadd.gs.washington.edu/download
cd $dbs
mkdir CADD
cd CADD
wget http://krishna.gs.washington.edu/download/CADD/v1.4/GRCh37/InDels.tsv.gz
wget http://krishna.gs.washington.edu/download/CADD/v1.4/GRCh37/whole_genome_SNVs.tsv.gz

#Install fathmm-MKL - https://github.com/HAShihab/fathmm-MKL
cd $dbs
mkdir fathmm-MKL
cd fathmm-MKL
wget http://fathmm.biocompute.org.uk/database/fathmm-MKL_Current.tab.gz
tabix -p bed fathmm-MKL_Current.tab.gz

#install OMIM (you might need a license - installation only possible after ngs-bits including NGSD is installed)
#cd $dbs
#mkdir OMIM
#cd OMIM
#manual download of ftp://ftp.omim.org/OMIM/genemap.txt
#manual download of ftp://ftp.omim.org/OMIM/mim2gene.txt
#php $src/Tools/db_converter_omim.php | $ngsbits/BedSort | bgzip > omim.bed.gz
#tabix -p bed omim.bed.gz

#Install HGMD (you need a license)
#manual download https://portal.biobase-international.com/cgi-bin/portal/login.cgi 
#cat HGMD_PRO_2018.2_hg19.vcf | php $src/Tools/db_converter_hgmd.php | bgzip > HGMD_PRO_2018_2_fixed.vcf.gz
#tabix -p vcf HGMD_PRO_2018_2_fixed.vcf.gz

