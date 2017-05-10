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

#install dbSNP
cd $dbs
mkdir dbSNP
cd dbSNP
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b149_GRCh37p13/VCF/00-All.vcf.gz
zcat 00-All.vcf.gz | php $src/Tools/db_converter_dbsnp.php | $vcflib/vcfbreakmulti | $ngsbits/VcfLeftNormalize -ref $genome | $ngsbits/VcfStreamSort | bgzip > dbsnp_b149.vcf.gz
tabix -p vcf dbsnp_b149.vcf.gz

#Install REPEATMASKER
cd $dbs
mkdir RepeatMasker
cd RepeatMasker
wget -O - http://www.repeatmasker.org/genomes/hg19/RepeatMasker-rm405-db20140131/hg19.fa.out.gz | gunzip > hg19.fa.out
perl $tools/RepeatMasker/util/rmOutToGFF3.pl hg19.fa.out > RepeatMasker.gff
cat RepeatMasker.gff | php $src/Tools/db_converter_repeatmasker.php | $ngsbits/BedSort > RepeatMasker.bed

#Install dbNSFP
cd $dbs
mkdir dbNSFP
cd dbNSFP
wget ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFPv2.9.2.zip
unzip dbNSFPv2.9.2.zip
head -n 1 dbNSFP2.9.2_variant.chr1 > dbNSFPv2.9.2.txt
cat dbNSFP2.9.2_variant.chr* | egrep -v "^#"  >> dbNSFPv2.9.2.txt
rm -rf dbNSFP2.9_gene.complete* dbNSFP2.9.2_variant* try* search*
bgzip dbNSFPv2.9.2.txt
tabix -s 1 -b 2 -e 2 dbNSFPv2.9.2.txt.gz

#Install 1000G
cd $dbs
mkdir 1000G
cd 1000G
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz
zcat ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz | $vcflib/vcfbreakmulti | $ngsbits/VcfLeftNormalize -ref $genome | $ngsbits/VcfStreamSort | bgzip > 1000g_v5b.vcf.gz
tabix -p vcf 1000g_v5b.vcf.gz

#Install ExAC
cd $dbs
mkdir ExAC
cd ExAC
wget ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz
zcat ExAC.r0.3.1.sites.vep.vcf.gz | $vcflib/vcfbreakmulti | $ngsbits/VcfLeftNormalize -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_exac.php | bgzip > ExAC_r0.3.1.vcf.gz
tabix -p vcf ExAC_r0.3.1.vcf.gz

#Install CLINVAR
cd $dbs
mkdir ClinVar
cd ClinVar
wget -O - ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/archive/2017/clinvar_20170130.vcf.gz | gunzip > clinvar_20170130.vcf
cat clinvar_20170130.vcf | php $src/Tools/db_converter_clinvar.php > clinvar_20170130_converted.vcf

#Install Kaviar
cd $dbs
mkdir Kaviar
cd Kaviar
wget http://s3-us-west-2.amazonaws.com/kaviar-160204-public/Kaviar-160204-Public-hg19-trim.vcf.tar
tar -xf Kaviar-160204-Public-hg19-trim.vcf.tar
mv Kaviar-160204-Public/vcfs/*.vcf.gz* .
rm -rf mv Kaviar-160204-Public *.vcf.tar
zcat Kaviar-160204-Public-hg19-trim.vcf.gz | $vcflib/vcfbreakmulti | $ngsbits/VcfLeftNormalize -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_kaviar.php | bgzip > Kaviar_160204.vcf.gz
tabix -p vcf Kaviar_160204.vcf.gz

#install OMIM (you might need a license)
#cd $dbs
#mkdir OMIM
#cd OMIM
#manual download ftp://ftp.omim.org/OMIM/genemap.txt
#manual download ftp://ftp.omim.org/OMIM/mim2gene.txt
#php $src/Tools/db_converter_omim.php > omim.bed

#Install HGMD (you need a license)
#manual download https://portal.biobase-international.com/cgi-bin/portal/login.cgi 
#cat HGMD_PRO_2016_4.vcf | php $src/Tools/db_converter_hgmd.php > HGMD_PRO_2016_4_fixed.vcf

#install COSMIC (you need a license)
#cd $dbs
#mkdir COSMIC
#cd COSMIC
#manual download http://cancer.sanger.ac.uk/files/cosmic/current_release/VCF/CosmicCodingMuts.vcf.gz
#manual download http://cancer.sanger.ac.uk/files/cosmic/current_release/VCF/CosmicNonCodingVariants.vcf.gz
#zcat CosmicCodingMuts.vcf.gz CosmicNonCodingVariants.vcf.gz | php $src/Tools/db_converter_cosmic.php | $ngsbits/VcfStreamSort -a > cosmic.vcf
#bgzip cosmic.vcf
#tabix -p vcf cosmic.vcf.gz