#!/bin/bash
set -e
set -o pipefail
set -o verbose

root=`pwd`
src=$root/../src/
tools=$root/tools/
dbs=$root/dbs/
ngsbits=$tools/ngs-bits/bin
genome=$root/genomes/GRCh37.fa

#Install CancerHotspots.org
cd $dbs
mkdir -p cancerhotspots
cd cancerhotspots
wget https://www.cancerhotspots.org/files/hotspots_v2.xls
wget http://download.cbioportal.org/cancerhotspots/cancerhotspots.v2.maf.gz
ssconvert -O 'separator="	" format=raw' -T Gnumeric_stf:stf_assistant -S hotspots_v2.xls hotspots.tsv
php $src/Tools/db_converter_cancerhotspots.php -in hotspots.tsv.0 -maf cancerhotspots.v2.maf.gz -out cancerhotspots_snv.tsv
rm hotspots_v2.xls
rm hotspots.tsv.0 
rm hotspots.tsv.1
rm cancerhotspots.v2.maf.gz

#Install ClinGen dosage sensitivity - ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/clingen
cd $dbs
mkdir ClinGen
cd ClinGen
wget ftp://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh37.tsv
cat ClinGen_gene_curation_list_GRCh37.tsv | php $src/Tools/db_converter_clingen_dosage.php > dosage_sensitive_disease_genes.bed
$ngsbits/BedSort -in dosage_sensitive_disease_genes.bed -out dosage_sensitive_disease_genes.bed

#Install NCG6.0 - information about oncogenes and tumor suppressor genes
cd $dbs
mkdir NCG6.0
cd NCG6.0
curl --silent --request POST --url http://ncg.kcl.ac.uk/download.php --data "filename=NCG6_tsgoncogene.tsv&downloadtsgoncogene=Download" --output NCG6.0_oncogene.tsv

#Install REPEATMASKER - http://www.repeatmasker.org/species/hg.html
cd $dbs
mkdir RepeatMasker
cd RepeatMasker
wget -O - http://www.repeatmasker.org/genomes/hg19/RepeatMasker-rm405-db20140131/hg19.fa.out.gz | gunzip > hg19.fa.out
cat hg19.fa.out | php $src/Tools/db_converter_repeatmasker.php | $ngsbits/BedSort | bgzip > RepeatMasker.bed.gz
tabix -p bed RepeatMasker.bed.gz

#Install ClinVar - https://www.ncbi.nlm.nih.gov/clinvar/
cd $dbs
mkdir ClinVar
cd ClinVar
wget -O - ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/archive_2.0/2021/clinvar_20210110.vcf.gz | gunzip | php $src/Tools/db_converter_clinvar.php | bgzip > clinvar_20210110_converted.vcf.gz
tabix -p vcf clinvar_20210110_converted.vcf.gz
#CNVs
wget -O - ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/variant_summary_2021-01.txt.gz | gunzip > variant_summary_2021-01.txt
cat variant_summary_2021-01.txt | php $src/Tools/db_converter_clinvar_cnvs.php 5 "Pathogenic/Likely pathogenic" | sort | uniq > clinvar_cnvs_2021-01.bed
$ngsbits/BedSort -with_name -in clinvar_cnvs_2021-01.bed -out clinvar_cnvs_2021-01.bed

#Install HGNC - ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/
cd $dbs
mkdir HGNC
cd HGNC
wget -O - ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt > hgnc_complete_set.tsv
wget -O - ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/withdrawn.txt > hgnc_withdrawn.tsv

#Install gnomAD (genome data) - http://gnomad.broadinstitute.org/downloads
cd $dbs
mkdir gnomAD
cd gnomAD
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.1.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -compression_level 0 -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php -header > gnomAD_genome_r2.1.1.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.2.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -compression_level 0 -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.3.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -compression_level 0 -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.4.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -compression_level 0 -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.5.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -compression_level 0 -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.6.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -compression_level 0 -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.7.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -compression_level 0 -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.8.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -compression_level 0 -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.9.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -compression_level 0 -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.10.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -compression_level 0 -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.11.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -compression_level 0 -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.12.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -compression_level 0 -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.13.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -compression_level 0 -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.14.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -compression_level 0 -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.15.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -compression_level 0 -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.16.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -compression_level 0 -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.17.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -compression_level 0 -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.18.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -compression_level 0 -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.19.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -compression_level 0 -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.20.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -compression_level 0 -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.21.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -compression_level 0 -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.22.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -compression_level 0 -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.X.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -compression_level 0 -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1.vcf
bgzip gnomAD_genome_r2.1.1.vcf
tabix -p vcf gnomAD_genome_r2.1.1.vcf.gz

#Install phyloP for VEP - https://www.ensembl.org/info/docs/tools/vep/script/vep_example.html#gerp
cd $dbs
mkdir phyloP
cd phyloP
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/phyloP100way/hg19.100way.phyloP100way.bw

#Install CADD for VEP - http://cadd.gs.washington.edu/download
cd $dbs
mkdir CADD
cd CADD
wget -O - http://krishna.gs.washington.edu/download/CADD/v1.6/GRCh37/InDels.tsv.gz > CADD_InDels_1.6.tsv.gz
wget -O - http://krishna.gs.washington.edu/download/CADD/v1.6/GRCh37/InDels.tsv.gz.tbi > CADD_InDels_1.6.tsv.gz.tbi
wget -O - http://krishna.gs.washington.edu/download/CADD/v1.6/GRCh37/whole_genome_SNVs.tsv.gz > CADD_SNVs_1.6.tsv.gz
wget -O - http://krishna.gs.washington.edu/download/CADD/v1.6/GRCh37/whole_genome_SNVs.tsv.gz.tbi > CADD_SNVs_1.6.tsv.gz.tbi
zcat CADD_InDels_1.6.tsv.gz | php $src/Tools/db_converter_cadd.php -in - -out - | VcfStreamSort | bgzip > CADD_InDels_1.6.vcf.gz
tabix -f -p vcf CADD_InDels_1.6.vcf.gz
zcat CADD_SNVs_1.6.tsv.gz | php $src/Tools/db_converter_cadd.php -in - -out - | VcfStreamSort | bgzip > CADD_SNVs_1.6.vcf.gz
tabix -f -p vcf CADD_SNVs_1.6.vcf.gz
VcfCheck -in CADD_InDels_1.6.vcf.gz -lines 0
VcfCheck -in CADD_SNVs_1.6.vcf.gz -lines 0

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

#GiaB NA12878 reference data
cd $dbs
mkdir -p GIAB/NA12878
cd GIAB/NA12878
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz -O high_conf_variants.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz.tbi -O high_conf_variants.vcf.gz.tbi
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed -O high_conf_regions.bed

#download annotation file for SpliceAI
cd $dbs
mkdir SpliceAI
cd SpliceAI
wget https://download.imgag.de/ahsturm1/spliceai_scores_2021_02_03.vcf.gz -O spliceai_scores_2021_02_03.vcf.gz
tabix -p vcf spliceai_scores_2021_02_03.vcf.gz

#download annotation file for MMSplice
cd $dbs
mkdir MMSplice
cd MMSplice
wget https://download.imgag.de/ahsturm1/mmsplice_scores_2021_02_03.vcf.gz -O mmsplice_scores_2021_02_03.vcf.gz
tabix -p vcf mmsplice_scores_2021_02_03.vcf.gz

#install OMIM (you might need a license, installation only possible after ngs-bits including NGSD is installed)
#cd $dbs
#mkdir OMIM
#cd OMIM
#manual download of ftp://ftp.omim.org/OMIM/genemap2.txt
#php $src/Tools/db_converter_omim.php | $ngsbits/BedSort -with_name > omim.bed
#bgzip -c omim.bed > omim.bed.gz
#tabix -p bed omim.bed.gz

#Install HGMD (you need a license)
#manual download of files hgmd_pro_2020.4_hg19.vcf and hgmd_pro-2020.4.dump.gz from https://portal.biobase-international.com/cgi-bin/portal/login.cgi 
#cat hgmd_pro_2020.4_hg19.vcf | php $src/Tools/db_converter_hgmd.php | bgzip > HGMD_PRO_2020_4_fixed.vcf.gz
#tabix -p vcf HGMD_PRO_2020_4_fixed.vcf.gz
##CNVs
#zcat hgmd_pro-2020.4.dump.gz | php $src/Tools/db_converter_hgmd_cnvs.php > HGMD_CNVS_2020_4.bed
#$ngsbits/BedSort -with_name -in HGMD_CNVS_2020_4.bed -out HGMD_CNVS_2020_4.bed


#Install COSMIC Cancer Mutation Census CMC  (you need a license, CMC tsv.gz file has to be downloaded manually from https://cancer.sanger.ac.uk/cmc/download)
#cd $dbs
#mkdir -p COSMIC
#cd COSMIC
##Login and download cmc_export.v92.tsv.gz from https://cancer.sanger.ac.uk/cmc/download to data/dbs/COSMIC. There is no download API for CMC file.
#mv cmc_export.v92.tsv.gz cmc_export.v92.tar.gz #CMC file is incorrectly name as tsv.gz when downloaded from COSMIC
#tar -xOzf cmc_export.v92.tar.gz cmc_export.tsv | php $src/Tools/db_converter_cosmic.php -in - -out cmc_export.vcf.gz

#install NGSD
#
#The usage of the NGSD annotation is optional.
#To generate the required VCF and BED files follow the instructions at https://github.com/imgag/ngs-bits/blob/master/doc/install_ngsd.md#export-ngsd-annotation-data (Export NGSD annotation data)
#The generated files have to be linked to "$data_folder/dbs/NGSD/" as symbolic links and have to be named as follows:
#	- "NGSD_germline.vcf.gz" for the germline export
#	- "NGSD_somatic.vcf.gz" for the somatic export
#	- "NGSD_genes.bed" for the gene info
#It is required the these files are symbolic links to avoid wrong annotations while performing a new export! (megSAP will check if these files are symlinks and fail if not)
#The actual files should be updated on regular bases (e.g. by using a cron job).
#Example code (generates a date based subfolder and links the generated files to the main folder):
#cd $dbs
#curdate=`date +"%Y-%m-%d"`
#mkdir $curdate
#cd $curdate
#$ngsbits/NGSDExportAnnotationData -variants NGSD_germline_unsorted.vcf -genes NGSD_genes.bed
#$ngsbits/VcfStreamSort -in NGSD_germline_unsorted.vcf -out NGSD_germline.vcf
#bgzip -c NGSD_germline.vcf > NGSD_germline.vcf.gz
#tabix -p vcf NGSD_germline.vcf.gz
#rm NGSD_germline_unsorted.vcf
#rm NGSD_germline.vcf
#$ngsbits/NGSDExportAnnotationData -variants NGSD_somatic_unsorted.vcf -mode somatic
#$ngsbits/VcfStreamSort -in NGSD_somatic_unsorted.vcf -out NGSD_somatic.vcf
#bgzip -c NGSD_somatic.vcf > NGSD_somatic.vcf.gz
#tabix -p vcf NGSD_somatic.vcf.gz
#rm NGSD_somatic_unsorted.vcf
#rm NGSD_somatic.vcf
#cd ..
#
#rm -f NGSD_germline.vcf.gz.tbi
#rm -f NGSD_somatic.vcf.gz.tbi
#rm -f NGSD_germline.vcf.gz
#rm -f NGSD_genes.bed
#rm -f NGSD_somatic.vcf.gz
#ln -s $curdate/NGSD_genes.bed NGSD_genes.bed
#ln -s $curdate/NGSD_germline.vcf.gz NGSD_germline.vcf.gz
#ln -s $curdate/NGSD_somatic.vcf.gz NGSD_somatic.vcf.gz
#ln -s $curdate/NGSD_germline.vcf.gz.tbi NGSD_germline.vcf.gz.tbi
#ln -s $curdate/NGSD_somatic.vcf.gz.tbi NGSD_somatic.vcf.gz.tbi


#Install fathmm-XF for AIdiva (custom VEP annotation) - https://fathmm.biocompute.org.uk/fathmm-xf
cd $dbs
mkdir fathmm-XF
cd fathmm-XF
wget http://fathmm.biocompute.org.uk/fathmm-xf/fathmm_xf_coding.vcf.gz
python3 $src/Tools/create_fathmm-XF_vcf.py fathmm_xf_coding.vcf.gz hg19_fathmm_xf_coding.vcf
bgzip hg19_fathmm_xf_coding.vcf
tabix -p vcf hg19_fathmm_xf_coding.vcf.gz
rm fathmm_xf_coding.vcf.gz

#Install ABB score for AIdiva (custom VEP annotation) - https://github.com/Francesc-Muyas/ABB
cd $dbs
mkdir ABB
cd ABB
python3 $src/Tools/create_ABB-SCORE_bed.py ABB_SCORE.txt hg19_ABB-SCORE.bed
bgzip hg19_ABB-SCORE.bed
tabix -p bed hg19_ABB-SCORE.bed.gz
rm ABB_SCORE.txt

#Install Eigen phred for AIdiva (custom VEP annotation)
cd $dbs
mkdir Eigen
cd Eigen
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr1.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr2.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr3.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr4.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr5.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr6.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr7.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr8.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr9.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr10.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr11.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr12.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr13.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr14.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr15.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr16.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr17.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr18.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr19.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr20.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr21.gz
wget -c http://web.corral.tacc.utexas.edu/WGSAdownload/resources/Eigen/Eigen_hg19_combined.tab.chr22.gz
python3 $src/Tools/create_eigen_vcf.py $(ls -m *.tab.chr*.gz | tr -d '[:space:]') hg19_Eigen-phred_coding_chrom1-22.vcf
#python3 $src/Tools/create_eigen_vcf-files.py $(ls -m *.tab.chr*.gz | tr -d '[:space:]') hg19_Eigen-phred_coding_chrom1-22.vcf hg19_Eigen-phred_noncoding_chrom1-22.vcf
bgzip hg19_Eigen-phred_coding_chrom1-22.vcf
#bgzip hg19_Eigen-phred_noncoding_chrom1-22.vcf
tabix -p vcf hg19_Eigen-phred_coding_chrom1-22.vcf.gz
#tabix -p vcf hg19_Eigen-phred_noncoding_chrom1-22.vcf.gz
rm *.tab.chr*.gz

#Install Condel for AIdiva (custom VEP annotation)
cd $dbs
mkdir Condel
cd Condel
wget -c https://bbglab.irbbarcelona.org/fannsdb/downloads/fannsdb.tsv.gz
python3 $src/Tools/create_Condel_vcf.py fannsdb.tsv.gz hg19_precomputed_Condel.vcf
bgzip hg19_precomputed_Condel.vcf
tabix -p vcf hg19_precomputed_Condel.vcf.gz
rm fannsdb.tsv.gz

cd $dbs
mkdir MutationAssessor
cd MutationAssessor
wget -c http://mutationassessor.org/r3/MA_scores_rel3_hg19_full.tar.gz
tar -xzf MA_scores_rel3_hg19_full.tar.gz
python3 $src/Tools/create_MutationAssessor_vcf.py $(ls -m MA_scores_rel3_hg19_full/*.csv | tr -d '[:space:]') hg19_precomputed_MutationAssessor.vcf
bgzip hg19_precomputed_MutationAssessor.vcf
tabix -p vcf hg19_precomputed_MutationAssessor.vcf.gz
rm MA_scores_rel3_hg19_full.tar.gz
rm -r MA_scores_rel3_hg19_full/

#Install segment duplication and simple repeats for AIdiva (custom VEP annotation)
cd $dbs
mkdir UCSC
cd UCSC
wget -O hg19_genomicSuperDups.txt.gz ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz
gunzip hg19_genomicSuperDups.txt.gz | cut -d$'\t' -f2,3,4,27 > hg19_genomicSuperDups.bed
grep -v "#" hg19_genomicSuperDups.bed | sort -k1,1 -k2,2n -k3,3n -t$'\t' | bgzip -c > hg19_genomicSuperDups.bed.gz
tabix -p bed hg19_genomicSuperDups.bed.gz
wget -O hg19_simpleRepeat.txt.gz ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz
gunzip hg19_simpleRepeat.txt.gz
cut -f2,3,4,11 hg19_simpleRepeat.txt > hg19_simpleRepeat.bed
grep -v '#' hg19_simpleRepeat.bed | sort -k1,1 -k2,2n -k3,3n -t '\t' | bgzip -c > hg19_simpleRepeat.bed.gz
tabix -p bed hg19_simpleRepeat.bed.gz

#Install phastCons46way vertebrate for AIdiva (custom VEP annotation)
cd $dbs
mkdir -p phastCons
cd phastCons
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr1.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr10.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr11.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr11_gl000202_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr12.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr13.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr14.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr15.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr16.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr17.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr17_ctg5_hap1.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr17_gl000203_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr17_gl000204_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr17_gl000205_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr17_gl000206_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr18.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr18_gl000207_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr19.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr19_gl000208_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr19_gl000209_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr1_gl000191_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr1_gl000192_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr2.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr20.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr21.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr21_gl000210_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr22.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr3.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr4.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr4_ctg9_hap1.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr4_gl000193_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr4_gl000194_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr5.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr6.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr6_apd_hap1.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr6_cox_hap2.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr6_dbb_hap3.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr6_mann_hap4.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr6_mcf_hap5.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr6_qbl_hap6.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr6_ssto_hap7.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr7.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr7_gl000195_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr8.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr8_gl000196_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr8_gl000197_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr9.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr9_gl000198_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr9_gl000199_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr9_gl000200_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chr9_gl000201_random.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrM.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000211.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000212.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000213.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000214.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000215.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000216.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000217.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000218.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000219.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000220.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000221.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000222.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000223.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000224.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000225.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000226.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000227.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000228.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000229.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000230.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000231.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000232.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000233.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000234.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000235.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000236.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000237.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000238.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000239.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000240.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000241.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000242.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000243.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000244.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000245.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000246.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000247.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000248.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrUn_gl000249.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrX.phastCons46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/vertebrate/chrY.phastCons46way.wigFix.gz
wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
wget -c http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
wget -c http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigCat
chmod +x wigToBigWig
chmod +x bigWigCat
gunzip *.wigFix.gz
for f in *.wigFix; do ./wigToBigWig -fixedSummaries -keepAllChromosomes $f hg19.chrom.sizes $(basename $f ".wigFix").bw; done;
rm *.wigFix
./bigWigCat hg19_phastCons46way_vertebrate.bw *.phastCons46way.bw
rm *.phastCons46way.bw

#Install phastCons46way primate for AIdiva (custom VEP annotation)
cd $dbs
mkdir -p phastCons
cd phastCons
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr1.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr10.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr11.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr11_gl000202_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr12.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr13.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr14.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr15.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr16.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr17.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr17_ctg5_hap1.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr17_gl000203_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr17_gl000204_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr17_gl000205_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr17_gl000206_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr18.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr18_gl000207_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr19.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr19_gl000208_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr19_gl000209_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr1_gl000191_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr1_gl000192_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr2.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr20.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr21.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr21_gl000210_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr22.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr3.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr4.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr4_ctg9_hap1.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr4_gl000193_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr4_gl000194_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr5.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr6.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr6_apd_hap1.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr6_cox_hap2.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr6_dbb_hap3.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr6_mann_hap4.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr6_mcf_hap5.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr6_qbl_hap6.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr6_ssto_hap7.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr7.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr7_gl000195_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr8.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr8_gl000196_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr8_gl000197_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr9.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr9_gl000198_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr9_gl000199_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr9_gl000200_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chr9_gl000201_random.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrM.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000211.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000212.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000213.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000214.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000215.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000216.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000217.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000218.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000219.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000220.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000221.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000222.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000223.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000224.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000225.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000226.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000227.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000228.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000229.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000230.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000231.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000232.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000233.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000234.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000235.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000236.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000237.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000238.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000239.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000240.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000241.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000242.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000243.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000244.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000245.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000246.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000247.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000248.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrUn_gl000249.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrX.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/chrY.phastCons46way.primates.wigFix.gz
wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
wget -c http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
wget -c http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigCat
chmod +x wigToBigWig
chmod +x bigWigCat
gunzip *.wigFix.gz
for f in *.wigFix; do ./wigToBigWig -fixedSummaries -keepAllChromosomes $f hg19.chrom.sizes $(basename $f ".wigFix").bw; done;
rm *.wigFix
./bigWigCat hg19_phastCons46way_primate.bw *.phastCons46way.primates.bw
rm *.phastCons46way.primates.bw

#Install phastCons46way mammal for AIdiva (custom VEP annotation)
cd $dbs
mkdir -p phastCons
cd phastCons
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr1.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr10.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr11.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr11_gl000202_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr12.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr13.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr14.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr15.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr16.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr17.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr17_ctg5_hap1.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr17_gl000203_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr17_gl000204_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr17_gl000205_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr17_gl000206_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr18.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr18_gl000207_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr19.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr19_gl000208_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr19_gl000209_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr1_gl000191_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr1_gl000192_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr2.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr20.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr21.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr21_gl000210_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr22.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr3.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr4.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr4_ctg9_hap1.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr4_gl000193_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr4_gl000194_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr5.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr6.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr6_apd_hap1.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr6_cox_hap2.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr6_dbb_hap3.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr6_mann_hap4.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr6_mcf_hap5.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr6_qbl_hap6.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr6_ssto_hap7.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr7.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr7_gl000195_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr8.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr8_gl000196_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr8_gl000197_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr9.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr9_gl000198_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr9_gl000199_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr9_gl000200_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chr9_gl000201_random.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrM.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000211.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000212.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000213.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000214.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000215.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000216.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000217.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000218.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000219.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000220.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000221.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000222.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000223.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000224.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000225.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000226.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000227.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000228.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000229.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000230.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000231.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000232.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000233.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000234.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000235.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000236.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000237.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000238.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000239.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000240.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000241.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000242.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000243.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000244.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000245.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000246.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000247.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000248.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrUn_gl000249.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrX.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/placentalMammals/chrY.phastCons46way.placental.wigFix.gz
wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
wget -c http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
wget -c http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigCat
chmod +x wigToBigWig
chmod +x bigWigCat
gunzip *.wigFix.gz
for f in *.wigFix; do ./wigToBigWig -fixedSummaries -keepAllChromosomes $$f hg19.chrom.sizes $$(basename $$f ".wigFix").bw; done;
rm *.wigFix
./bigWigCat hg19_phastCons46way_mammal.bw *.phastCons46way.placental.bw
rm *.phastCons46way.placental.bw

#Install phyloP46way primate for AIdiva (custom VEP annotation)
cd $dbs
mkdir -p phyloP
cd phyloP
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr1.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr10.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr11.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr11_gl000202_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr12.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr13.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr14.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr15.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr16.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr17.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr17_ctg5_hap1.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr17_gl000203_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr17_gl000204_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr17_gl000205_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr17_gl000206_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr18.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr18_gl000207_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr19.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr19_gl000208_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr19_gl000209_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr1_gl000191_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr1_gl000192_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr2.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr20.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr21.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr21_gl000210_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr22.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr3.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr4.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr4_ctg9_hap1.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr4_gl000193_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr4_gl000194_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr5.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr6.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr6_apd_hap1.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr6_cox_hap2.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr6_dbb_hap3.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr6_mann_hap4.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr6_mcf_hap5.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr6_qbl_hap6.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr6_ssto_hap7.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr7.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr7_gl000195_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr8.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr8_gl000196_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr8_gl000197_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr9.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr9_gl000198_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr9_gl000199_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr9_gl000200_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chr9_gl000201_random.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrM.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000211.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000212.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000213.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000214.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000215.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000216.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000217.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000218.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000219.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000220.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000221.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000222.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000223.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000224.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000225.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000227.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000228.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000229.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000230.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000231.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000232.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000233.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000234.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000235.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000236.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000237.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000238.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000239.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000240.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000241.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000242.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000243.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000244.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000245.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000246.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000247.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000248.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrUn_gl000249.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrX.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/primates/chrY.phyloP46way.primate.wigFix.gz
wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
wget -c http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
wget -c http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigCat
chmod +x wigToBigWig
chmod +x bigWigCat
gunzip *.wigFix.gz
for f in *.wigFix; do ./wigToBigWig -fixedSummaries -keepAllChromosomes $f hg19.chrom.sizes $(basename $f ".wigFix").bw; done;
rm *.wigFix
./bigWigCat hg19_phyloP46way_primate.bw *.phyloP46way.primate.bw
rm *.phyloP46way.primate.bw

#Install phyloP46way mammal for AIdiva (custom VEP annotation)
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr1.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr10.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr11.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr11_gl000202_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr12.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr13.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr14.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr15.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr16.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr17.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr17_ctg5_hap1.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr17_gl000203_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr17_gl000204_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr17_gl000205_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr17_gl000206_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr18.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr18_gl000207_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr19.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr19_gl000208_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr19_gl000209_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr1_gl000191_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr1_gl000192_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr2.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr20.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr21.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr21_gl000210_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr22.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr3.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr4.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr4_ctg9_hap1.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr4_gl000193_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr4_gl000194_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr5.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr6.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr6_apd_hap1.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr6_cox_hap2.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr6_dbb_hap3.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr6_mann_hap4.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr6_mcf_hap5.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr6_qbl_hap6.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr6_ssto_hap7.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr7.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr7_gl000195_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr8.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr8_gl000196_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr8_gl000197_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr9.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr9_gl000198_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr9_gl000199_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr9_gl000200_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chr9_gl000201_random.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrM.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000211.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000212.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000213.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000214.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000215.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000216.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000217.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000218.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000219.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000220.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000221.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000222.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000223.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000224.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000225.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000227.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000228.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000229.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000230.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000231.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000232.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000233.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000234.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000235.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000236.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000237.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000238.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000239.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000240.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000241.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000242.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000243.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000244.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000245.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000246.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000247.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000248.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrUn_gl000249.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrX.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/placentalMammals/chrY.phyloP46way.placental.wigFix.gz
wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
wget -c http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
wget -c http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigCat
chmod +x wigToBigWig
chmod +x bigWigCat
gunzip *.wigFix.gz
for f in *.wigFix; do ./wigToBigWig -fixedSummaries -keepAllChromosomes $f hg19.chrom.sizes $(basename $f ".wigFix").bw; done;
rm *.wigFix
./bigWigCat hg19_phyloP46way_mammal.bw *.phyloP46way.placental.bw
rm *.phyloP46way.placental.bw

#Install phyloP46way vertebrate for AIdiva (custom VEP annotation)
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr1.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr10.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr11.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr11_gl000202_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr12.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr13.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr14.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr15.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr16.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr17.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr17_ctg5_hap1.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr17_gl000203_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr17_gl000204_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr17_gl000205_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr17_gl000206_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr18.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr18_gl000207_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr19.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr19_gl000208_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr19_gl000209_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr1_gl000191_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr1_gl000192_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr2.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr20.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr21.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr21_gl000210_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr22.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr3.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr4.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr4_ctg9_hap1.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr4_gl000193_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr4_gl000194_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr5.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr6.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr6_apd_hap1.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr6_cox_hap2.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr6_dbb_hap3.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr6_mann_hap4.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr6_mcf_hap5.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr6_qbl_hap6.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr6_ssto_hap7.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr7.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr7_gl000195_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr8.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr8_gl000196_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr8_gl000197_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr9.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr9_gl000198_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr9_gl000199_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr9_gl000200_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chr9_gl000201_random.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrM.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000211.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000212.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000213.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000214.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000215.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000216.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000217.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000218.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000219.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000220.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000221.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000222.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000223.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000224.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000225.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000227.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000228.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000229.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000230.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000231.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000232.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000233.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000234.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000235.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000236.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000237.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000238.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000239.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000240.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000241.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000242.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000243.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000244.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000245.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000246.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000247.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000248.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrUn_gl000249.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrX.phyloP46way.wigFix.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP46way/vertebrate/chrY.phyloP46way.wigFix.gz
wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
wget -c http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
wget -c http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigCat
chmod +x wigToBigWig
chmod +x bigWigCat
gunzip *.wigFix.gz
for f in *.wigFix; do ./wigToBigWig -fixedSummaries -keepAllChromosomes $f hg19.chrom.sizes $(basename $f ".wigFix").bw; done;
rm *.wigFix
./bigWigCat hg19_phyloP46way_vertebrate.bw *.phyloP46way.bw
rm *.phyloP46way.bw
