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

#Download ensembl transcripts database
cd $dbs
mkdir -p Ensembl
cd Ensembl
wget https://ftp.ensembl.org/pub/release-110/gff3/homo_sapiens/Homo_sapiens.GRCh38.110.gff3.gz
gunzip Homo_sapiens.GRCh38.110.gff3.gz

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

#Install ClinGen dosage sensitivity - http://ftp.ncbi.nlm.nih.gov/pub/dbVar/clingen
cd $dbs
mkdir ClinGen
cd ClinGen
wget http://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv
cat ClinGen_gene_curation_list_GRCh38.tsv | php $src/Tools/db_converter_clingen_dosage.php > dosage_sensitive_disease_genes_GRCh38.bed
$ngsbits/BedSort -in dosage_sensitive_disease_genes_GRCh38.bed -out dosage_sensitive_disease_genes_GRCh38.bed

#Install NCG7.1 - information about oncogenes and tumor suppressor genes
cd $dbs
mkdir NCG7.1
cd NCG7.1
curl 'http://network-cancer-genes.org/download.php' --silent -X POST --data-raw 'downloadcancergenes=Download' | sed -e "1s/^/#/" > ncg.tsv
php $src/Tools/db_converter_ncg.php -in ncg.tsv -prefix NCG7.1 -outfolder "."

#Install REPEATMASKER - http://www.repeatmasker.org/species/hg.html
cd $dbs
mkdir RepeatMasker
cd RepeatMasker
wget -O - http://www.repeatmasker.org/genomes/hg38/RepeatMasker-rm405-db20140131/hg38.fa.out.gz | gunzip > hg38.fa.out
cat hg38.fa.out | php $src/Tools/db_converter_repeatmasker.php | $ngsbits/BedSort > RepeatMasker_GRCh38.bed

#Install ClinVar - https://www.ncbi.nlm.nih.gov/clinvar/
cd $dbs
mkdir ClinVar 
cd ClinVar
wget -O - http://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2024/clinvar_20240127.vcf.gz | gunzip | php $src/Tools/db_converter_clinvar.php | $ngsbits/VcfStreamSort | bgzip > clinvar_20240127_converted_GRCh38.vcf.gz
tabix -C -m 9 -p vcf clinvar_20240127_converted_GRCh38.vcf.gz
#CNVs
wget -O - http://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/variant_summary_2023-11.txt.gz | gunzip > variant_summary_2023-11.txt
cat variant_summary_2023-11.txt | php $src/Tools/db_converter_clinvar_cnvs.php 5 "Pathogenic/Likely pathogenic" | sort | uniq > clinvar_cnvs_2024-02.bed
$ngsbits/BedSort -with_name -in clinvar_cnvs_2024-02.bed -out clinvar_cnvs_2024-02.bed

#Install HGNC - http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/
cd $dbs
mkdir HGNC
cd HGNC
wget -O - http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt > hgnc_complete_set.tsv
wget -O - http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/withdrawn.txt > hgnc_withdrawn.tsv

#Install gnomAD (genome data) - 
cd $dbs
mkdir gnomAD
cd gnomAD
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr1.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php -header > gnomAD_genome_v3.1.2_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr2.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_v3.1.2_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr3.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_v3.1.2_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr4.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_v3.1.2_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr5.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_v3.1.2_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr6.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_v3.1.2_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr7.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_v3.1.2_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr8.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_v3.1.2_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr9.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_v3.1.2_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr10.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_v3.1.2_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr11.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_v3.1.2_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr12.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_v3.1.2_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr13.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_v3.1.2_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr14.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_v3.1.2_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr15.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_v3.1.2_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr16.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_v3.1.2_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr17.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_v3.1.2_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr18.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_v3.1.2_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr19.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_v3.1.2_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr20.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_v3.1.2_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr21.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_v3.1.2_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr22.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_v3.1.2_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chrX.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_v3.1.2_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chrY.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | php $src/Tools/db_converter_gnomad.php >> gnomAD_genome_v3.1.2_GRCh38.vcf
bgzip gnomAD_genome_v3.1.2_GRCh38.vcf
tabix -C -m 9 -p vcf gnomAD_genome_v3.1.2_GRCh38.vcf.gz
wget -O - https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1/vcf/genomes/gnomad.genomes.v3.1.sites.chrM.vcf.bgz | gunzip | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort >> gnomAD_genome_v3.1.mito_GRCh38.vcf
bgzip gnomAD_genome_v3.1.mito_GRCh38.vcf
tabix -C -m 9 -p vcf gnomAD_genome_v3.1.mito_GRCh38.vcf.gz

#Install phyloP
cd $dbs
mkdir phyloP
cd phyloP
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw

#Install CADD
cd $dbs
mkdir CADD
cd CADD
wget -O - http://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.indel.tsv.gz > CADD_InDels_1.6_GRCh38.tsv.gz
wget -O - http://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz > CADD_SNVs_1.6_GRCh38.tsv.gz
zcat CADD_InDels_1.6_GRCh38.tsv.gz | php $src/Tools/db_converter_cadd.php -build GRCh38 -in - -out - | $ngsbits/VcfStreamSort | bgzip > CADD_InDels_1.6_GRCh38.vcf.gz
tabix -f -C -m 9 -p vcf CADD_InDels_1.6_GRCh38.vcf.gz
zcat CADD_SNVs_1.6_GRCh38.tsv.gz | php $src/Tools/db_converter_cadd.php -build GRCh38 -in - -out - | $ngsbits/VcfStreamSort | bgzip > CADD_SNVs_1.6_GRCh38.vcf.gz
tabix -f -C -m 9 -p vcf CADD_SNVs_1.6_GRCh38.vcf.gz
$ngsbits/VcfCheck -in CADD_InDels_1.6_GRCh38.vcf.gz -lines 1000 -ref $genome
$ngsbits/VcfCheck -in CADD_SNVs_1.6_GRCh38.vcf.gz -lines 1000 -ref $genome

#download and convert REVEL - https://sites.google.com/site/revelgenomics/downloads
cd $dbs
mkdir REVEL
cd REVEL
wget https://zenodo.org/record/7072866/files/revel-v1.3_all_chromosomes.zip
unzip -p revel-v1.3_all_chromosomes.zip | php $src/Tools/db_converter_revel.php > tmp.vcf
$ngsbits/VcfSort -in tmp.vcf -out REVEL_1.3.vcf
rm tmp.vcf
bgzip REVEL_1.3.vcf
tabix -f -C -m 9 -p vcf REVEL_1.3.vcf.gz
$ngsbits/VcfCheck -in REVEL_1.3.vcf.gz -lines 1000 -ref $genome

#download and convert AlphaMissense - Attention: for non-commercial use only!
cd $dbs
mkdir AlphaMissense
cd AlphaMissense
wget https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz
php $src/Tools/db_converter_alphamissense.php AlphaMissense_hg38.tsv.gz > AlphaMissense_hg38.vcf
$ngsbits/VcfSort -in AlphaMissense_hg38.vcf -out AlphaMissense_hg38.vcf
bgzip AlphaMissense_hg38.vcf
tabix -p vcf AlphaMissense_hg38.vcf.gz

#download annotation file for SpliceAI
cd $dbs
mkdir SpliceAI
cd SpliceAI
wget https://download.imgag.de/public/splicing/spliceai_scores_2023_12_20_GRCh38.vcf.gz -O spliceai_scores_2023_12_20_GRCh38.vcf.gz
tabix -C -m 9 -p vcf spliceai_scores_2023_12_20_GRCh38.vcf.gz

#download sniffles 
cd $dbs
mkdir -p TandemRepeats
cd TandemRepeats
wget https://github.com/PacificBiosciences/pbsv/raw/master/annotations/human_GRCh38_no_alt_analysis_set.trf.bed


#install OMIM (you might need a license, only possible after ngs-bits is installed - including reference genome and NGSD setup)
#cd $dbs
#mkdir OMIM
#cd OMIM
#manual download of http://ftp.omim.org/OMIM/genemap2.txt
#php $src/Tools/db_converter_omim.php | $ngsbits/BedSort -with_name > omim.bed

#Install HGMD (you need a license, only possible after ngs-bits is installed - including reference genome and NGSD setup)
#manual download of files HGMD_Pro_2023.3_hg38.vcf.gz  and hgmd_pro-2023.2.dump.gz from https://apps.ingenuity.com/ingsso/login
#zcat HGMD_Pro_2023.3_hg38.vcf.gz | php $src/Tools/db_converter_hgmd.php | bgzip > HGMD_PRO_2023_3_fixed.vcf.gz
#tabix -p vcf HGMD_PRO_2023_3_fixed.vcf.gz
##CNVs
#zcat hgmd_pro-2023.3.dump.gz | php $src/Tools/db_converter_hgmd_cnvs.php > HGMD_CNVS_2023_3.bed
#$ngsbits/BedSort -with_name -in HGMD_CNVS_2023_3.bed -out HGMD_CNVS_2023_3.bed

#Install COSMIC Cancer Mutation Census CMC  (you need a license, the files have to be downloaded manually from https://apps.ingenuity.com/ingsso/login)
#the necessary files are: CancerMutationCensus_AllData_Tsv_v98_GRCh38.tar, Cosmic_GenomeScreensMutant_Vcf_v98_GRCh38.tar, Cosmic_CompleteTargetedScreensMutant_Vcf_v98_GRCh38.tar, Cosmic_NonCodingVariants_Vcf_v98_GRCh38.tar
# unpack the files to get the necessary vcf.gz files: 
#cd $dbs
#mkdir -p COSMIC
#cd COSMIC
##Login and download the "Cancer Mutation Census Data": CMC.tar file from https://cancer.sanger.ac.uk/cmc/download. Extract the file and move the cmc_export.tsv.gz to data/dbs/COSMIC. There is no download API for CMC file.
#gunzip -c CancerMutationCensus_AllData_v98_GRCh38.tsv.gz | php db_converter_cosmic.php -in_cmc cmc_export.tsv -in_genome_vcf Cosmic_GenomeScreensMutant_v98_GRCh38.vcf.gz -in_non_coding_vcf Cosmic_NonCodingVariants_v98_GRCh38.vcf.gz -in_target_screens_vcf Cosmic_CompleteTargetedScreensMutant_v98_GRCh38.vcf.gz -out cmc_export_v98.vcf.gz
#install NGSD
#
#The usage of the NGSD annotation is optional. 
#To generate the required VCF, BEDPE and BED files follow the instructions at https://github.com/imgag/ngs-bits/blob/master/doc/install_ngsd.md#export-ngsd-annotation-data (Export NGSD annotation data)
#The generated files have to be linked to "$data_folder/dbs/NGSD/" as symbolic links and have to be named as follows:
#	- "NGSD_germline.vcf.gz" for the germline export 
#	- "NGSD_somatic.vcf.gz" for the somatic export 
#	- "NGSD_genes.bed" for the gene info
#	- "sv_deletion.bedpe.gz" for deletions
#	- "sv_duplication.bedpe.gz" for duplications
#	- "sv_insertion.bedpe.gz" for insertions
#	- "sv_inversion.bedpe.gz" for inversions
#	- "sv_translocation.bedpe.gz" for translocation
#	- "sv_breakpoint_density.igv" for breakpoints
#It is required the these files are symbolic links to avoid wrong annotations while performing a new export (megSAP will check if these files are symlinks and fail if not)
#The actual files should be updated on regular bases (e.g. by using a cron job).
#Example code (generates a date based subfolder and links the generated files to the main folder):
#cd $dbs
#curdate=`date +"%Y-%m-%d"`
#mkdir $curdate
#cd $curdate
#$ngsbits/NGSDExportAnnotationData -germline NGSD_germline_unsorted.vcf -somatic NGSD_somatic_unsorted.vcf -genes NGSD_genes.bed
#$ngsbits/VcfStreamSort -in NGSD_germline_unsorted.vcf -out NGSD_germline.vcf
#bgzip -c NGSD_germline.vcf > NGSD_germline.vcf.gz
#tabix -p vcf NGSD_germline.vcf.gz
#rm NGSD_germline_unsorted.vcf
#rm NGSD_germline.vcf
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
#
#bgzip -c $curdate/sv_deletion.bedpe > $curdate/sv_deletion.bedpe.gz
#bgzip -c $curdate/sv_duplication.bedpe > $curdate/sv_duplication.bedpe.gz
#bgzip -c $curdate/sv_insertion.bedpe > $curdate/sv_insertion.bedpe.gz
#bgzip -c $curdate/sv_inversion.bedpe > $curdate/sv_inversion.bedpe.gz
#$ngsbits/BedpeSort -in $curdate/sv_translocation.bedpe -out $curdate/sv_translocation.bedpe
#bgzip -c $curdate/sv_translocation.bedpe > $curdate/sv_translocation.bedpe.gz
#
#tabix -0 -b 2 -e 5 $curdate/sv_deletion.bedpe.gz
#tabix -0 -b 2 -e 5 $curdate/sv_duplication.bedpe.gz
#tabix -0 -b 2 -e 3 $curdate/sv_insertion.bedpe.gz
#tabix -0 -b 2 -e 5 $curdate/sv_inversion.bedpe.gz
#tabix -0 -b 2 -e 3 $curdate/sv_translocation.bedpe.gz
#
#rm $curdate/sv_*.bedpe
#
#rm -f sv_deletion.bedpe.gz
#rm -f sv_duplication.bedpe.gz
#rm -f sv_insertion.bedpe.gz
#rm -f sv_inversion.bedpe.gz
#rm -f sv_translocation.bedpe.gz
#rm -f sv_breakpoint_density.igv
#
#rm -f sv_deletion.bedpe.gz.tbi
#rm -f sv_duplication.bedpe.gz.tbi
#rm -f sv_insertion.bedpe.gz.tbi
#rm -f sv_inversion.bedpe.gz.tbi
#rm -f sv_translocation.bedpe.gz.tbi
#
#ln -s $curdate/sv_deletion.bedpe.gz sv_deletion.bedpe.gz
#ln -s $curdate/sv_duplication.bedpe.gz sv_duplication.bedpe.gz
#ln -s $curdate/sv_insertion.bedpe.gz sv_insertion.bedpe.gz
#ln -s $curdate/sv_inversion.bedpe.gz sv_inversion.bedpe.gz
#ln -s $curdate/sv_translocation.bedpe.gz sv_translocation.bedpe.gz
#
#ln -s $curdate/sv_deletion.bedpe.gz.tbi sv_deletion.bedpe.gz.tbi
#ln -s $curdate/sv_duplication.bedpe.gz.tbi sv_duplication.bedpe.gz.tbi
#ln -s $curdate/sv_insertion.bedpe.gz.tbi sv_insertion.bedpe.gz.tbi
#ln -s $curdate/sv_inversion.bedpe.gz.tbi sv_inversion.bedpe.gz.tbi
#ln -s $curdate/sv_translocation.bedpe.gz.tbi sv_translocation.bedpe.gz.tbi
#
#ln -s $curdate/sv_breakpoint_density.igv sv_breakpoint_density.igv
