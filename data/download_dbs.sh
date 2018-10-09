#!/bin/bash
set -e
set -o pipefail
set -o verbose

root=`pwd`
src=$root/../src/
tools=$root/tools/
dbs=$root/dbs/
ngsbits=$tools/ngs-bits/bin
vcflib=$tools/vcflib/bin
genome=$root/genomes/GRCh37.fa

#Install REPEATMASKER - http://www.repeatmasker.org/species/hg.html
cd $dbs
mkdir RepeatMasker
cd RepeatMasker
wget -O - http://www.repeatmasker.org/genomes/hg19/RepeatMasker-rm405-db20140131/hg19.fa.out.gz | gunzip > hg19.fa.out
perl $tools/RepeatMasker/util/rmOutToGFF3.pl hg19.fa.out > RepeatMasker.gff
cat RepeatMasker.gff | php $src/Tools/db_converter_repeatmasker.php | $ngsbits/BedSort | bgzip > RepeatMasker.bed
tabix -p bed RepeatMasker.bed.gz
rm -rf hg19.fa.out RepeatMasker.gff

#Install ClinVar - https://www.ncbi.nlm.nih.gov/clinvar/
cd $dbs
mkdir ClinVar
cd ClinVar
wget -O - ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/archive_2.0/2018/clinvar_20180805.vcf.gz | gunzip | php $src/Tools/db_converter_clinvar.php | bgzip > clinvar_20180805_converted.vcf.gz
tabix -p vcf clinvar_20180805_converted.vcf.gz

#Install gnomAD (genome data) - http://gnomad.broadinstitute.org/downloads
cd $dbs
cd gnomAD
wget -O - https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr1.vcf.bgz | gunzip | $vcflib/vcfbreakmulti | $ngsbits/VcfLeftNormalize -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php -header > gnomAD_genome_r2.0.2.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr2.vcf.bgz | gunzip | $vcflib/vcfbreakmulti | $ngsbits/VcfLeftNormalize -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.0.2.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr3.vcf.bgz | gunzip | $vcflib/vcfbreakmulti | $ngsbits/VcfLeftNormalize -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.0.2.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr4.vcf.bgz | gunzip | $vcflib/vcfbreakmulti | $ngsbits/VcfLeftNormalize -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.0.2.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr5.vcf.bgz | gunzip | $vcflib/vcfbreakmulti | $ngsbits/VcfLeftNormalize -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.0.2.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr6.vcf.bgz | gunzip | $vcflib/vcfbreakmulti | $ngsbits/VcfLeftNormalize -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.0.2.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr7.vcf.bgz | gunzip | $vcflib/vcfbreakmulti | $ngsbits/VcfLeftNormalize -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.0.2.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr8.vcf.bgz | gunzip | $vcflib/vcfbreakmulti | $ngsbits/VcfLeftNormalize -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.0.2.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr9.vcf.bgz | gunzip | $vcflib/vcfbreakmulti | $ngsbits/VcfLeftNormalize -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.0.2.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr10.vcf.bgz | gunzip | $vcflib/vcfbreakmulti | $ngsbits/VcfLeftNormalize -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.0.2.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr11.vcf.bgz | gunzip | $vcflib/vcfbreakmulti | $ngsbits/VcfLeftNormalize -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.0.2.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr12.vcf.bgz | gunzip | $vcflib/vcfbreakmulti | $ngsbits/VcfLeftNormalize -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.0.2.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr13.vcf.bgz | gunzip | $vcflib/vcfbreakmulti | $ngsbits/VcfLeftNormalize -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.0.2.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr14.vcf.bgz | gunzip | $vcflib/vcfbreakmulti | $ngsbits/VcfLeftNormalize -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.0.2.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr15.vcf.bgz | gunzip | $vcflib/vcfbreakmulti | $ngsbits/VcfLeftNormalize -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.0.2.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr16.vcf.bgz | gunzip | $vcflib/vcfbreakmulti | $ngsbits/VcfLeftNormalize -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.0.2.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr17.vcf.bgz | gunzip | $vcflib/vcfbreakmulti | $ngsbits/VcfLeftNormalize -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.0.2.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr18.vcf.bgz | gunzip | $vcflib/vcfbreakmulti | $ngsbits/VcfLeftNormalize -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.0.2.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr19.vcf.bgz | gunzip | $vcflib/vcfbreakmulti | $ngsbits/VcfLeftNormalize -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.0.2.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr20.vcf.bgz | gunzip | $vcflib/vcfbreakmulti | $ngsbits/VcfLeftNormalize -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.0.2.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr21.vcf.bgz | gunzip | $vcflib/vcfbreakmulti | $ngsbits/VcfLeftNormalize -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.0.2.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr22.vcf.bgz | gunzip | $vcflib/vcfbreakmulti | $ngsbits/VcfLeftNormalize -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.0.2.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chrX.vcf.bgz | gunzip | $vcflib/vcfbreakmulti | $ngsbits/VcfLeftNormalize -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.0.2.vcf
bgzip gnomAD_genome_r2.0.2.vcf
tabix -p vcf gnomAD_genome_r2.0.2.vcf.gz

#Install phyloP for VEP - https://www.ensembl.org/info/docs/tools/vep/script/vep_example.html#gerp
cd $dbs
mkdir phyloP
cd phyloP
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/phyloP100way/hg19.100way.phyloP100way.bw

#Install CADD for VEP - http://cadd.gs.washington.edu/download
cd $dbs
mkdir CADD
cd CADD
wget http://krishna.gs.washington.edu/download/CADD/v1.4/GRCh37/InDels.tsv.gz
wget http://krishna.gs.washington.edu/download/CADD/v1.4/GRCh37/whole_genome_SNVs.tsv.gz

#Install fathmm-MKL for VEP - https://github.com/HAShihab/fathmm-MKL
cd $dbs
mkdir fathmm-MKL
cd fathmm-MKL
wget http://fathmm.biocompute.org.uk/database/fathmm-MKL_Current.tab.gz
tabix -p bed fathmm-MKL_Current.tab.gz

#Install REVEL for VEP - https://sites.google.com/site/revelgenomics/downloads
cd $dbs
mkdir REVEL
cd REVEL
wget https://rothsj06.u.hpc.mssm.edu/revel/revel_all_chromosomes.csv.zip
unzip -p revel_all_chromosomes.csv.zip | tr ',' '\t' | sed '1s/.*/#&/' | bgzip > revel_all_chromosomes.tsv.gz
tabix -f -s 1 -b 2 -e 2 revel_all_chromosomes.tsv.gz

#Install dbscSNV for VEP - https://academic.oup.com/nar/article/42/22/13534/2411339
cd $dbs
mkdir dbscSNV
cd dbscSNV
wget ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbscSNV1.1.zip
unzip dbscSNV1.1.zip
head -n1 dbscSNV1.1.chr1 > h
cat dbscSNV1.1.chr* | grep -v ^chr | cat h - | bgzip -c > dbscSNV1.1_GRCh37.txt.gz
tabix -s 1 -b 2 -e 2 -c c dbscSNV1.1_GRCh37.txt.gz

#install OMIM (you might need a license - installation only possible after ngs-bits including NGSD is installed)
#cd $dbs
#mkdir OMIM
#cd OMIM
#manual download of ftp://ftp.omim.org/OMIM/genemap.txt
#manual download of ftp://ftp.omim.org/OMIM/mim2gene.txt
#php $src/Tools/db_converter_omim.php | $ngsbits/BedSort > omim.bed
#bgzip -c omim.bed > omim.bed.gz
#tabix -p bed omim.bed.gz

#Install HGMD (you need a license)
#manual download https://portal.biobase-international.com/cgi-bin/portal/login.cgi 
#cat HGMD_PRO_2018.3_hg19.vcf | php $src/Tools/db_converter_hgmd.php | bgzip > HGMD_PRO_2018_3_fixed.vcf.gz
#tabix -p vcf HGMD_PRO_2018_3_fixed.vcf.gz

