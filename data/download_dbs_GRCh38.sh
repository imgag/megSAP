#!/bin/bash
set -e
set -o pipefail
set -o verbose

root=`pwd`
src=$root/../src/
tools=$root/tools/
dbs=$root/dbs/
ngsbits=$tools/ngs-bits/bin
genome=$root/genomes/GRCh38.fa

# TODO: Find GRCh38 version or create liftover
# #Install CancerHotspots.org
# cd $dbs
# mkdir -p cancerhotspots
# cd cancerhotspots
# wget https://www.cancerhotspots.org/files/hotspots_v2.xls
# wget http://download.cbioportal.org/cancerhotspots/cancerhotspots.v2.maf.gz
# ssconvert -O 'separator="	" format=raw' -T Gnumeric_stf:stf_assistant -S hotspots_v2.xls hotspots.tsv
# php $src/Tools/db_converter_cancerhotspots.php -in hotspots.tsv.0 -maf cancerhotspots.v2.maf.gz -out cancerhotspots_snv.tsv
# rm hotspots_v2.xls
# rm hotspots.tsv.0 
# rm hotspots.tsv.1
# rm cancerhotspots.v2.maf.gz

#Install ClinGen dosage sensitivity - ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/clingen
cd $dbs
mkdir ClinGen
cd ClinGen
wget ftp://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv
cat ClinGen_gene_curation_list_GRCh38.tsv | php $src/Tools/db_converter_clingen_dosage.php > dosage_sensitive_disease_genes_GRCh38.bed
$ngsbits/BedSort -in dosage_sensitive_disease_genes_GRCh38.bed -out dosage_sensitive_disease_genes_GRCh38.bed

# TODO: Find GRCh38 version or create liftover
# #Install NCG6.0 - information about oncogenes and tumor suppressor genes
# cd $dbs
# mkdir NCG6.0
# cd NCG6.0
# curl --silent --request POST --url http://ncg.kcl.ac.uk/download.php --data "filename=NCG6_tsgoncogene.tsv&downloadtsgoncogene=Download" --output NCG6.0_oncogene.tsv

#Install REPEATMASKER - http://www.repeatmasker.org/species/hg.html
cd $dbs
mkdir RepeatMasker
cd RepeatMasker
wget -O - http://www.repeatmasker.org/genomes/hg38/RepeatMasker-rm405-db20140131/hg38.fa.out.gz | gunzip > hg38.fa.out
cat hg38.fa.out | php $src/Tools/db_converter_repeatmasker.php | $ngsbits/BedSort | bgzip > RepeatMasker_GRCh38.bed.gz
tabix -p bed RepeatMasker_GRCh38.bed.gz

#Install ClinVar - https://www.ncbi.nlm.nih.gov/clinvar/
cd $dbs
mkdir ClinVar
cd ClinVar
wget -O - ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2021/clinvar_20210424.vcf.gz | gunzip | php $src/Tools/db_converter_clinvar.php | bgzip > clinvar_20210424_converted_GRCh38.vcf.gz
tabix -p vcf clinvar_20210424_converted_GRCh38.vcf.gz
#CNVs
wget -O - ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/variant_summary_2021-04.txt.gz | gunzip > variant_summary_2021-04.txt
cat variant_summary_2021-04.txt | php $src/Tools/db_converter_clinvar_cnvs.php 5 "Pathogenic/Likely pathogenic" | sort | uniq > clinvar_cnvs_2021-04.bed
$ngsbits/BedSort -with_name -in clinvar_cnvs_2021-04.bed -out clinvar_cnvs_2021-04.bed

#Install HGNC - ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/
cd $dbs
mkdir HGNC
cd HGNC
wget -O - ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt > hgnc_complete_set.tsv
wget -O - ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/withdrawn.txt > hgnc_withdrawn.tsv

#Install gnomAD (genome data) - 
cd $dbs
mkdir gnomAD
cd gnomAD
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.1.liftover_grch38.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php -header > gnomAD_genome_r2.1.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.2.liftover_grch38.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.3.liftover_grch38.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.4.liftover_grch38.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.5.liftover_grch38.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.6.liftover_grch38.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.7.liftover_grch38.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.8.liftover_grch38.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.9.liftover_grch38.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.10.liftover_grch38.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.11.liftover_grch38.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.12.liftover_grch38.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.13.liftover_grch38.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.14.liftover_grch38.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.15.liftover_grch38.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.16.liftover_grch38.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.17.liftover_grch38.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.18.liftover_grch38.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.19.liftover_grch38.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.20.liftover_grch38.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.21.liftover_grch38.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.22.liftover_grch38.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gnomad-public/release/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.X.liftover_grch38.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_r2.1.1_GRCh38.vcf
bgzip gnomAD_genome_r2.1.1_GRCh38.vcf
tabix -p vcf gnomAD_genome_r2.1.1_GRCh38.vcf.gz

#Install phyloP for VEP - https://www.ensembl.org/info/docs/tools/vep/script/vep_example.html#gerp
cd $dbs
mkdir phyloP
cd phyloP
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw

#Install CADD for VEP - http://cadd.gs.washington.edu/download
cd $dbs
mkdir CADD
cd CADD
wget -O - http://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.indel.tsv.gz > CADD_InDels_1.6_GRCh38.tsv.gz
wget -O - http://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.indel.tsv.gz.tbi > CADD_InDels_1.6_GRCh38.tsv.gz.tbi
wget -O - http://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz > CADD_SNVs_1.6_GRCh38.tsv.gz
wget -O - http://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz.tbi > CADD_SNVs_1.6_GRCh38.tsv.gz.tbi
zcat CADD_InDels_1.6_GRCh38.tsv.gz | php $src/Tools/db_converter_cadd.php -build GRCh38 -in - -out - | $ngsbits/VcfStreamSort | bgzip > CADD_InDels_1.6_GRCh38.vcf.gz
tabix -f -p vcf CADD_InDels_1.6_GRCh38.vcf.gz
zcat CADD_SNVs_1.6_GRCh38.tsv.gz | php $src/Tools/db_converter_cadd.php -build GRCh38 -in - -out - | $ngsbits/VcfStreamSort | bgzip > CADD_SNVs_1.6_GRCh38.vcf.gz
tabix -f -p vcf CADD_SNVs_1.6_GRCh38.vcf.gz
$ngsbits/VcfCheck -in CADD_InDels_1.6_GRCh38.vcf.gz -lines 0 –ref $genome
$ngsbits/VcfCheck -in CADD_SNVs_1.6_GRCh38.vcf.gz -lines 0 –ref $genome

# TODO: Find GRCh38 version or create liftover
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
wget https://rothsj06.u.hpc.mssm.edu/revel_grch38_all_chromosomes.csv.zip
unzip -p revel_grch38_all_chromosomes.csv.zip | tr ',' '\t' | sed '1s/.*/#&/' | bgzip > revel_grch38_all_chromosomes.tsv.gz
tabix -f -s 1 -b 2 -e 2 revel_grch38_all_chromosomes.tsv.gz

#Install dbscSNV for VEP - https://academic.oup.com/nar/article/42/22/13534/2411339
cd $dbs
mkdir dbscSNV
cd dbscSNV
wget ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbscSNV1.1.zip
unzip dbscSNV1.1.zip
head -n1 dbscSNV1.1.chr1 > h
# cat dbscSNV1.1.chr* | grep -v ^chr | cat h - | bgzip -c > dbscSNV1.1_GRCh37.txt.gz
# tabix -s 1 -b 2 -e 2 -c c dbscSNV1.1_GRCh37.txt.gz
cat dbscSNV1.1.chr* | grep -v ^chr | sort -k5,5 -k6,6n | cat h - | awk '$5 != "."' | bgzip -c > dbscSNV1.1_GRCh38.txt.gz
tabix -s 5 -b 6 -e 6 -c c dbscSNV1.1_GRCh38.txt.gz



#GiaB NA12878 reference data
cd $dbs
mkdir -p GIAB/NA12878
cd GIAB/NA12878
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz -O high_conf_variants_GRCh38.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz.tbi -O high_conf_variants_GRCh38.vcf.gz.tbi
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed -O high_conf_regions_GRCh38.bed

# TODO: Create GRCh38 version
#download annotation file for SpliceAI
# cd $dbs
# mkdir SpliceAI
# cd SpliceAI
# wget https://download.imgag.de/ahsturm1/spliceai_scores_2021_02_03.vcf.gz -O spliceai_scores_2021_02_03.vcf.gz
# tabix -p vcf spliceai_scores_2021_02_03.vcf.gz

# TODO: Create GRCh38 version
#download annotation file for MMSplice
# cd $dbs
# mkdir MMSplice
# cd MMSplice
# wget https://download.imgag.de/ahsturm1/mmsplice_scores_2021_02_03.vcf.gz -O mmsplice_scores_2021_02_03.vcf.gz
# tabix -p vcf mmsplice_scores_2021_02_03.vcf.gz

#install OMIM (you might need a license, only possible after ngs-bits is installed - including reference genome and NGSD setup)
#cd $dbs
#mkdir OMIM
#cd OMIM
#manual download of ftp://ftp.omim.org/OMIM/genemap2.txt
#php $src/Tools/db_converter_omim.php | $ngsbits/BedSort -with_name > omim.bed
#bgzip -c omim.bed > omim.bed.gz
#tabix -p bed omim.bed.gz

#Install HGMD (you need a license, only possible after ngs-bits is installed - including reference genome and NGSD setup)
#manual download of files hgmd_pro_2021.1_hg19.vcf and hgmd_pro-2021.1.dump.gz from https://apps.ingenuity.com/ingsso/login
#cat hgmd_pro_2021.1_hg19.vcf | php $src/Tools/db_converter_hgmd.php | bgzip > HGMD_PRO_2021_1_fixed.vcf.gz
#tabix -p vcf HGMD_PRO_2021_1_fixed.vcf.gz
##CNVs
#zcat hgmd_pro-2020.4.dump.gz | php $src/Tools/db_converter_hgmd_cnvs.php > HGMD_CNVS_2021_1.bed
#$ngsbits/BedSort -with_name -in HGMD_CNVS_2021_1.bed -out HGMD_CNVS_2021_1.bed


#Install COSMIC Cancer Mutation Census CMC  (you need a license, CMC tsv.gz file has to be downloaded manually from https://cancer.sanger.ac.uk/cmc/download)
#cd $dbs
#mkdir -p COSMIC
#cd COSMIC
##Login and download cmc_export.v92.tsv.gz from https://cancer.sanger.ac.uk/cmc/download to data/dbs/COSMIC. There is no download API for CMC file.
#mv cmc_export.v92.tsv.gz cmc_export.v92.tar.gz #CMC file is incorrectly name as tsv.gz when downloaded from COSMIC
#tar -xOzf cmc_export.v92.tar.gz cmc_export.tsv | php $src/Tools/db_converter_cosmic.php -build GRCh38 -in - -out cmc_export.vcf.gz

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
