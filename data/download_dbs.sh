#!/bin/bash
set -e
set -o pipefail
set -o verbose

root=`pwd`
src=$root/../src/

# Get data_folder path
SETTINGS_FILE=$root/../settings.ini
if [ ! -f "$SETTINGS_FILE" ]; then
    SETTINGS_FILE="$root/../settings.ini.default"
fi
DATA_FOLDER=$(grep -E "^data_folder" "$SETTINGS_FILE" | awk -F ' = ' '{print $2}' | sed "s|\[path\]|$(dirname "$root")|")

dbs=$DATA_FOLDER/dbs/
misc=$DATA_FOLDER/misc/
genome_dir=$DATA_FOLDER/genomes/
genome=$genome_dir/GRCh38.fa

mkdir -p $dbs
mkdir -p $misc

# Get ngs-bits container path
SETTINGS_FILE=$root/../settings.ini
if [ ! -f "$SETTINGS_FILE" ]; then
    SETTINGS_FILE="$root/../settings.ini.default"
fi
CONTAINER_FOLDER=$(grep -E "^container_folder" "$SETTINGS_FILE" | awk -F ' = ' '{print $2}' | sed "s|\[path\]|$(dirname "$root")|")
NGSBITS_VERSION=$(grep -E "^container_ngs-bits" "$SETTINGS_FILE" | awk -F ' = ' '{print $2}')
HTSLIB_VERSION=$(grep -E "^container_htslib" "$SETTINGS_FILE" | awk -F ' = ' '{print $2}')
MSISENSOR_VERSION=$(grep -E "^container_msisensor-pro" "$SETTINGS_FILE" | awk -F ' = ' '{print $2}')
htslib=$CONTAINER_FOLDER/htslib_$HTSLIB_VERSION.sif
ngsbits=$CONTAINER_FOLDER/ngs-bits_$NGSBITS_VERSION.sif
msisensor=$CONTAINER_FOLDER/msisensor-pro_$MSISENSOR_VERSION.sif

#Download ensembl transcripts database
cd $dbs
mkdir -p Ensembl
cd Ensembl
wget -O Homo_sapiens.GRCh38.112.gff3.gz https://ftp.ensembl.org/pub/release-112/gff3/homo_sapiens/Homo_sapiens.GRCh38.112.gff3.gz
gunzip Homo_sapiens.GRCh38.112.gff3.gz
wget -O Homo_sapiens.GRCh38.112.gtf.gz https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz
gunzip Homo_sapiens.GRCh38.112.gtf.gz
# create sorted & indexed file for methylartist
(grep ^"#" Homo_sapiens.GRCh38.112.gtf; grep -v ^"#" Homo_sapiens.GRCh38.112.gtf | sort -k1,1 -k4,4n | sed -e 's/^/chr/') | apptainer exec $htslib bgzip  > Homo_sapiens.GRCh38.112.gtf.gz
apptainer exec $htslib tabix -p gff Homo_sapiens.GRCh38.112.gtf.gz

#Download RefSeq transcripts database
cd $dbs
mkdir -p RefSeq
cd RefSeq
wget -O GCF_000001405.40_GRCh38.p14_genomic.gff.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz
zcat GCF_000001405.40_GRCh38.p14_genomic.gff.gz > Homo_sapiens.GRCh38.p14.gff3

#Install CancerHotspots.org
cd $dbs
mkdir -p cancerhotspots
cd cancerhotspots
wget -O hotspots_v2.xls https://www.cancerhotspots.org/files/hotspots_v2.xls
wget -O cancerhotspots.v2.maf.gz https://cbioportal-download.s3.amazonaws.com/cancerhotspots.v2.maf.gz
ssconvert -O 'separator="	" format=raw' -T Gnumeric_stf:stf_assistant -S hotspots_v2.xls hotspots.tsv
php $src/Install/db_converter_cancerhotspots.php -in hotspots.tsv.0 -maf cancerhotspots.v2.maf.gz -out cancerhotspots_snv.tsv
rm hotspots_v2.xls
rm hotspots.tsv.0 
rm hotspots.tsv.1
rm cancerhotspots.v2.maf.gz

#Install ClinGen dosage sensitivity - http://ftp.ncbi.nlm.nih.gov/pub/dbVar/clingen
cd $dbs
mkdir -p ClinGen
cd ClinGen
wget -O ClinGen_gene_curation_list_GRCh38.tsv http://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv
cat ClinGen_gene_curation_list_GRCh38.tsv | php $src/Install/db_converter_clingen_dosage.php > dosage_sensitive_disease_genes_GRCh38.bed
apptainer exec $ngsbits BedSort -in dosage_sensitive_disease_genes_GRCh38.bed -out dosage_sensitive_disease_genes_GRCh38.bed

#Install NCG7.1 - information about oncogenes and tumor suppressor genes
cd $dbs
mkdir -p NCG7.1
cd NCG7.1
curl 'http://network-cancer-genes.org/download.php' --silent -X POST --data-raw 'downloadcancergenes=Download' | sed -e "1s/^/#/" > ncg.tsv
php $src/Install/db_converter_ncg.php -in ncg.tsv -prefix NCG7.1 -outfolder "."

#Install REPEATMASKER - http://www.repeatmasker.org/species/hg.html
cd $dbs
mkdir -p RepeatMasker
cd RepeatMasker
wget -O - http://www.repeatmasker.org/genomes/hg38/RepeatMasker-rm405-db20140131/hg38.fa.out.gz | gunzip > hg38.fa.out
cat hg38.fa.out | php $src/Install/db_converter_repeatmasker.php | apptainer exec $ngsbits BedSort > RepeatMasker_GRCh38.bed

#Install ClinVar - https://www.ncbi.nlm.nih.gov/clinvar/
cd $dbs
mkdir -p ClinVar 
cd ClinVar
wget -O - http://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2025/clinvar_20250128.vcf.gz | gunzip | php $src/Install/db_converter_clinvar.php | apptainer exec $ngsbits VcfStreamSort | apptainer exec $htslib bgzip > clinvar_20250128_converted_GRCh38.vcf.gz
apptainer exec $htslib tabix -C -m 9 -p vcf clinvar_20250128_converted_GRCh38.vcf.gz
#CNVs
wget -O - http://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/variant_summary_2025-02.txt.gz | gunzip > variant_summary_2025-02.txt
cat variant_summary_2025-02.txt | php $src/Install/db_converter_clinvar_cnvs.php 5 "Pathogenic/Likely pathogenic" | sort | uniq > clinvar_cnvs_2025-02.bed
apptainer exec $ngsbits BedSort -with_name -in clinvar_cnvs_2025-02.bed -out clinvar_cnvs_2025-02.bed

#Install HGNC - http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/
cd $dbs
mkdir -p HGNC
cd HGNC
wget -O - https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt > hgnc_complete_set.tsv
wget -O - https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/withdrawn.txt > hgnc_withdrawn.tsv

#Install gnomAD (genome data) - 
cd $dbs
mkdir -p gnomAD
cd gnomAD
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr1.vcf.bgz | gunzip | apptainer exec -B $genome_dir $ngsbits VcfLeftNormalize -stream -ref $genome | apptainer exec $ngsbits VcfStreamSort | php $src/Install/db_converter_gnomad.php -header > gnomAD_genome_v4.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr2.vcf.bgz | gunzip | apptainer exec -B $genome_dir $ngsbits VcfLeftNormalize -stream -ref $genome | apptainer exec $ngsbits VcfStreamSort | php $src/Install/db_converter_gnomad.php >> gnomAD_genome_v4.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr3.vcf.bgz | gunzip | apptainer exec -B $genome_dir $ngsbits VcfLeftNormalize -stream -ref $genome | apptainer exec $ngsbits VcfStreamSort | php $src/Install/db_converter_gnomad.php >> gnomAD_genome_v4.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr4.vcf.bgz | gunzip | apptainer exec -B $genome_dir $ngsbits VcfLeftNormalize -stream -ref $genome | apptainer exec $ngsbits VcfStreamSort | php $src/Install/db_converter_gnomad.php >> gnomAD_genome_v4.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr5.vcf.bgz | gunzip | apptainer exec -B $genome_dir $ngsbits VcfLeftNormalize -stream -ref $genome | apptainer exec $ngsbits VcfStreamSort | php $src/Install/db_converter_gnomad.php >> gnomAD_genome_v4.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr6.vcf.bgz | gunzip | apptainer exec -B $genome_dir $ngsbits VcfLeftNormalize -stream -ref $genome | apptainer exec $ngsbits VcfStreamSort | php $src/Install/db_converter_gnomad.php >> gnomAD_genome_v4.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr7.vcf.bgz | gunzip | apptainer exec -B $genome_dir $ngsbits VcfLeftNormalize -stream -ref $genome | apptainer exec $ngsbits VcfStreamSort | php $src/Install/db_converter_gnomad.php >> gnomAD_genome_v4.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr8.vcf.bgz | gunzip | apptainer exec -B $genome_dir $ngsbits VcfLeftNormalize -stream -ref $genome | apptainer exec $ngsbits VcfStreamSort | php $src/Install/db_converter_gnomad.php >> gnomAD_genome_v4.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr9.vcf.bgz | gunzip | apptainer exec -B $genome_dir $ngsbits VcfLeftNormalize -stream -ref $genome | apptainer exec $ngsbits VcfStreamSort | php $src/Install/db_converter_gnomad.php >> gnomAD_genome_v4.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr10.vcf.bgz | gunzip | apptainer exec -B $genome_dir $ngsbits VcfLeftNormalize -stream -ref $genome | apptainer exec $ngsbits VcfStreamSort | php $src/Install/db_converter_gnomad.php >> gnomAD_genome_v4.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr11.vcf.bgz | gunzip | apptainer exec -B $genome_dir $ngsbits VcfLeftNormalize -stream -ref $genome | apptainer exec $ngsbits VcfStreamSort | php $src/Install/db_converter_gnomad.php >> gnomAD_genome_v4.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr12.vcf.bgz | gunzip | apptainer exec -B $genome_dir $ngsbits VcfLeftNormalize -stream -ref $genome | apptainer exec $ngsbits VcfStreamSort | php $src/Install/db_converter_gnomad.php >> gnomAD_genome_v4.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr13.vcf.bgz | gunzip | apptainer exec -B $genome_dir $ngsbits VcfLeftNormalize -stream -ref $genome | apptainer exec $ngsbits VcfStreamSort | php $src/Install/db_converter_gnomad.php >> gnomAD_genome_v4.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr14.vcf.bgz | gunzip | apptainer exec -B $genome_dir $ngsbits VcfLeftNormalize -stream -ref $genome | apptainer exec $ngsbits VcfStreamSort | php $src/Install/db_converter_gnomad.php >> gnomAD_genome_v4.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr15.vcf.bgz | gunzip | apptainer exec -B $genome_dir $ngsbits VcfLeftNormalize -stream -ref $genome | apptainer exec $ngsbits VcfStreamSort | php $src/Install/db_converter_gnomad.php >> gnomAD_genome_v4.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr16.vcf.bgz | gunzip | apptainer exec -B $genome_dir $ngsbits VcfLeftNormalize -stream -ref $genome | apptainer exec $ngsbits VcfStreamSort | php $src/Install/db_converter_gnomad.php >> gnomAD_genome_v4.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr17.vcf.bgz | gunzip | apptainer exec -B $genome_dir $ngsbits VcfLeftNormalize -stream -ref $genome | apptainer exec $ngsbits VcfStreamSort | php $src/Install/db_converter_gnomad.php >> gnomAD_genome_v4.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr18.vcf.bgz | gunzip | apptainer exec -B $genome_dir $ngsbits VcfLeftNormalize -stream -ref $genome | apptainer exec $ngsbits VcfStreamSort | php $src/Install/db_converter_gnomad.php >> gnomAD_genome_v4.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr19.vcf.bgz | gunzip | apptainer exec -B $genome_dir $ngsbits VcfLeftNormalize -stream -ref $genome | apptainer exec $ngsbits VcfStreamSort | php $src/Install/db_converter_gnomad.php >> gnomAD_genome_v4.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr20.vcf.bgz | gunzip | apptainer exec -B $genome_dir $ngsbits VcfLeftNormalize -stream -ref $genome | apptainer exec $ngsbits VcfStreamSort | php $src/Install/db_converter_gnomad.php >> gnomAD_genome_v4.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr21.vcf.bgz | gunzip | apptainer exec -B $genome_dir $ngsbits VcfLeftNormalize -stream -ref $genome | apptainer exec $ngsbits VcfStreamSort | php $src/Install/db_converter_gnomad.php >> gnomAD_genome_v4.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr22.vcf.bgz | gunzip | apptainer exec -B $genome_dir $ngsbits VcfLeftNormalize -stream -ref $genome | apptainer exec $ngsbits VcfStreamSort | php $src/Install/db_converter_gnomad.php >> gnomAD_genome_v4.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chrX.vcf.bgz | gunzip | apptainer exec -B $genome_dir $ngsbits VcfLeftNormalize -stream -ref $genome | apptainer exec $ngsbits VcfStreamSort | php $src/Install/db_converter_gnomad.php >> gnomAD_genome_v4.1_GRCh38.vcf
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chrY.vcf.bgz | gunzip | apptainer exec -B $genome_dir $ngsbits VcfLeftNormalize -stream -ref $genome | apptainer exec $ngsbits VcfStreamSort | php $src/Install/db_converter_gnomad.php >> gnomAD_genome_v4.1_GRCh38.vcf
apptainer exec $htslib bgzip gnomAD_genome_v4.1_GRCh38.vcf
apptainer exec $htslib tabix -C -m 9 -p vcf gnomAD_genome_v4.1_GRCh38.vcf.gz
wget -O - https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/vcf/genomes/gnomad.genomes.v3.1.sites.chrM.vcf.bgz | gunzip | apptainer exec -B $genome_dir $ngsbits VcfLeftNormalize -stream -ref $genome | apptainer exec $ngsbits VcfStreamSort >> gnomAD_genome_v3.1.mito_GRCh38.vcf
apptainer exec $htslib bgzip gnomAD_genome_v3.1.mito_GRCh38.vcf
apptainer exec $htslib tabix -C -m 9 -p vcf gnomAD_genome_v3.1.mito_GRCh38.vcf.gz

#Install phyloP
cd $dbs
mkdir -p phyloP
cd phyloP
wget -O hg38.phyloP100way.bw http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw

#Install CADD
cd $dbs
mkdir -p CADD
cd CADD
wget -O - https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz > CADD_InDels_1.7_GRCh38.tsv.gz
wget -O - https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz > CADD_SNVs_1.7_GRCh38.tsv.gz
zcat CADD_InDels_1.7_GRCh38.tsv.gz | php $src/Install/db_converter_cadd.php -build GRCh38 | apptainer exec $ngsbits VcfStreamSort | apptainer exec $htslib bgzip > CADD_InDels_1.7_GRCh38.vcf.gz
apptainer exec $htslib tabix -f -C -m 9 -p vcf CADD_InDels_1.7_GRCh38.vcf.gz
zcat CADD_SNVs_1.7_GRCh38.tsv.gz | php $src/Install/db_converter_cadd.php -build GRCh38 | apptainer exec $ngsbits VcfStreamSort | apptainer exec $htslib bgzip > CADD_SNVs_1.7_GRCh38.vcf.gz
apptainer exec $htslib tabix -f -C -m 9 -p vcf CADD_SNVs_1.7_GRCh38.vcf.gz
apptainer exec -B $genome_dir $ngsbits VcfCheck -in CADD_InDels_1.7_GRCh38.vcf.gz -lines 1000 -ref $genome
apptainer exec -B $genome_dir $ngsbits VcfCheck -in CADD_SNVs_1.7_GRCh38.vcf.gz -lines 1000 -ref $genome
rm -rf CADD_SNVs_1.7_GRCh38.tsv.gz CADD_InDels_1.7_GRCh38.tsv.gz

#download and convert REVEL - https://sites.google.com/site/revelgenomics/downloads
cd $dbs
mkdir -p REVEL
cd REVEL
wget -O revel-v1.3_all_chromosomes.zip https://zenodo.org/record/7072866/files/revel-v1.3_all_chromosomes.zip
unzip -p revel-v1.3_all_chromosomes.zip | php $src/Install/db_converter_revel.php > tmp.vcf
apptainer exec $ngsbits VcfSort -in tmp.vcf -out REVEL_1.3.vcf
rm tmp.vcf
apptainer exec $htslib bgzip REVEL_1.3.vcf
apptainer exec $htslib tabix -f -C -m 9 -p vcf REVEL_1.3.vcf.gz
apptainer exec -B $genome_dir $ngsbits VcfCheck -in REVEL_1.3.vcf.gz -lines 1000 -ref $genome

#download and convert AlphaMissense - Attention: for non-commercial use only!
cd $dbs
mkdir -p AlphaMissense
cd AlphaMissense
wget -O AlphaMissense_hg38.tsv.gz https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz
php $src/Install/db_converter_alphamissense.php AlphaMissense_hg38.tsv.gz > AlphaMissense_hg38.vcf
apptainer exec $ngsbits VcfSort -in AlphaMissense_hg38.vcf -out AlphaMissense_hg38.vcf
apptainer exec $htslib bgzip AlphaMissense_hg38.vcf
apptainer exec $htslib tabix -p vcf AlphaMissense_hg38.vcf.gz

#download annotation file for SpliceAI
cd $dbs
mkdir -p SpliceAI
cd SpliceAI
wget https://megsap.de/download/SpliceAI/spliceai_scores_2024_08_26_GRCh38.vcf.gz -O spliceai_scores_2024_08_26_GRCh38.vcf.gz
apptainer exec $htslib tabix -C -m 9 -p vcf spliceai_scores_2024_08_26_GRCh38.vcf.gz

#download reference data for gene expression
cd $dbs
mkdir -p gene_expression
cd gene_expression
#change version number on update
wget -O - https://www.proteinatlas.org/download/tsv/rna_tissue_consensus.tsv.zip | gunzip > rna_tissue_consensus_v23.tsv

#download Ensembl data in GTF format - KEEP AT ENSEMBL VERSION 107, DB import for RNA works on Transcript base and will break if the transcripts change.
cd $dbs
mkdir -p gene_annotations
cd gene_annotations
wget -O - 'http://ftp.ensembl.org/pub/release-107/gtf/homo_sapiens/Homo_sapiens.GRCh38.107.gtf.gz' | gzip -cd | awk '{ if ($$1 !~ /^#/) { print "chr"$0 } else { print $0 } }' > GRCh38.gtf

#create hemoglobin FASTA file
cd $misc
wget -O - 'https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz' | zcat | \
	awk -v RS=">" -v FS="\n" '$$1 ~ / gene_symbol:(HBA1|HBA2|HBB) / { print ">"$$1; {for (i=2; i<=NF; i++) printf("%s", $$i)}; printf("\n") }' | \
	sed '/^>/s/ /|kraken:taxid|9606 /' \
	> human_hemoglobin_tx.fa

#download and normalize HG001/NA12878 reference data
mkdir -p $dbs/GIAB/NA12878/
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -O $dbs/GIAB/NA12878/high_conf_variants.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.bed -O $dbs/GIAB/NA12878/high_conf_regions.bed
zcat $dbs/GIAB/NA12878/high_conf_variants.vcf.gz | apptainer exec $ngsbits VcfBreakMulti | apptainer exec -B $genome_dir $ngsbits VcfFilter -remove_invalid -ref $genome | apptainer exec -B $genome_dir $ngsbits VcfLeftNormalize -stream -ref $genome | apptainer exec $ngsbits VcfStreamSort | apptainer exec $htslib bgzip > $dbs/GIAB/NA12878/high_conf_variants_normalized.vcf.gz
apptainer exec $htslib tabix $dbs/GIAB/NA12878/high_conf_variants_normalized.vcf.gz

#download and normalize HG002/NA24385 reference data
mkdir -p $dbs/GIAB/NA24385/
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -O $dbs/GIAB/NA24385/high_conf_variants.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -O $dbs/GIAB/NA24385/high_conf_regions.bed
zcat $dbs/GIAB/NA24385/high_conf_variants.vcf.gz | apptainer exec $ngsbits VcfBreakMulti | apptainer exec -B $genome_dir $ngsbits VcfFilter -remove_invalid -ref $genome | apptainer exec -B $genome_dir $ngsbits VcfLeftNormalize -stream -ref $genome | apptainer exec $ngsbits VcfStreamSort | apptainer exec $htslib bgzip > $dbs/GIAB/NA24385/high_conf_variants_normalized.vcf.gz
apptainer exec $htslib tabix $dbs/GIAB/NA24385/high_conf_variants_normalized.vcf.gz

#download and normalize HG002/NA24385 CMRG reference data
mkdir -p $dbs/GIAB/NA24385_CMRG/
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/SmallVariant/HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz -O $dbs/GIAB/NA24385_CMRG/high_conf_variants.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/SmallVariant/HG002_GRCh38_CMRG_smallvar_v1.00.bed -O $dbs/GIAB/NA24385_CMRG/high_conf_regions.bed
zcat $dbs/GIAB/NA24385_CMRG/high_conf_variants.vcf.gz | apptainer exec $ngsbits VcfBreakMulti | apptainer exec -B $genome_dir $ngsbits VcfFilter -remove_invalid -ref $genome | apptainer exec -B $genome_dir $ngsbits VcfLeftNormalize -stream -ref $genome | apptainer exec $ngsbits VcfStreamSort | apptainer exec $htslib bgzip > $dbs/GIAB/NA24385_CMRG/high_conf_variants_normalized.vcf.gz
apptainer exec $htslib tabix $dbs/GIAB/NA24385_CMRG/high_conf_variants_normalized.vcf.gz

#download reference genome for orad
cd $dbs
mkdir -p oradata
cd oradata
wget -O orad.2.6.1.tar.gz https://webdata.illumina.com/downloads/software/dragen-decompression/orad.2.6.1.tar.gz
tar xzf orad.2.6.1.tar.gz
rm orad.2.6.1.tar.gz
mv orad_2_6_1/oradata/refbin .
rm -rf orad_2_6_1

#create reference file for msisensor-pro
cd $dbs
mkdir -p msisensor-pro
cd msisensor-pro
apptainer exec -B $genome $msisensor msisensor-pro scan -d $genome -o msisensor_references_GRCh38.site

#download tandem-repeat file for sniffles2
cd $dbs
mkdir -p tandem-repeats
cd tandem-repeats
wget -O human_GRCh38_no_alt_analysis_set.trf.bed https://raw.githubusercontent.com/fritzsedlazeck/Sniffles/fdf6e6d334353a06872fe98f74fe68cc9a9a7d1f/annotations/human_GRCh38_no_alt_analysis_set.trf.bed

# # install OMIM (you might need a license; production NGSD has to be available and initialized)
# cd $dbs
# mkdir -p OMIM
# cd OMIM
# # manual download of http://ftp.omim.org/OMIM/genemap2.txt
# php $src/Install/db_converter_omim.php | apptainer exec $ngsbits BedSort -with_name > omim.bed
 
# # when using the containerized megSAP version:
# cd <path-to-host-data-folder>/dbs/
# mkdir -p OMIM
# cd OMIM
# # manual download of http://ftp.omim.org/OMIM/genemap2.txt
# apptainer exec -B <path-to-host-data-folder>:/megSAP/data/data_folder/,<path-to-settings.ini-with-NGSD-credentials>:/megSAP/settings.ini --pwd /megSAP/data/data_folder/dbs/OMIM megSAP_[version].sif php megSAP/src/Install/db_converter_omim.php | apptainer exec <downloaded-ngs-bits-container> BedSort -with_name > omim.bed


# # install HGMD (you need a license; production NGSD has to be available and initialized)
# cd $dbs
# mkdir -p HGMD
# cd HGMD
# # manual download of files HGMD_Pro_2024.4_hg38.vcf.gz and hgmd_pro-2024.4.dump.gz from https://apps.ingenuity.com/ingsso/login
# zcat HGMD_Pro_2024.4_hg38.vcf.gz | php $src/Install/db_converter_hgmd.php | apptainer exec $htslib bgzip > HGMD_PRO_2024_4_fixed.vcf.gz
# apptainer exec $htslib tabix -p vcf HGMD_PRO_2024_4_fixed.vcf.gz
# #CNVs
# zcat hgmd_pro-2024.4.dump.gz | php $src/Install/db_converter_hgmd_cnvs.php > HGMD_CNVS_2024_4.bed
# apptainer exec $ngsbits BedSort -with_name -in HGMD_CNVS_2024_4.bed -out HGMD_CNVS_2024_4.bed

# # when using the containerized megSAP version:
# cd <path-to-host-data-folder>/dbs/
# mkdir -p HGMD
# cd HGMD
# # manual download of files HGMD_Pro_2024.4_hg38.vcf.gz and hgmd_pro-2024.4.dump.gz from https://apps.ingenuity.com/ingsso/login
# apptainer exec -B <path-to-host-data-folder>:/megSAP/data/data_folder/ --pwd /megSAP/data/data_folder/dbs/HGMD megSAP_[version].sif sh -c "zcat HGMD_Pro_2024.4_hg38.vcf.gz | php /megSAP/src/Install/db_converter_hgmd.php | apptainer exec $htslib bgzip > HGMD_PRO_2024_4_fixed.vcf.gz"
# apptainer exec $htslib tabix -p vcf HGMD_PRO_2024_4_fixed.vcf.gz
# #CNVs
# apptainer exec -B <path-to-host-data-folder>:/megSAP/data/data_folder/,<path-to-settings.ini-with-NGSD-credentials>:/megSAP/settings.ini --pwd /megSAP/data/data_folder/dbs/HGMD megSAP_[version].sif sh -c "zcat hgmd_pro-2024.4.dump.gz | php /megSAP/src/Install/db_converter_hgmd_cnvs.php > HGMD_CNVS_2024_4.bed"
# apptainer exec <downloaded-ngs-bits-container> BedSort -with_name -in HGMD_CNVS_2024_4.bed -out HGMD_CNVS_2024_4.bed


# # install COSMIC Cancer Mutation Census CMC (you need a license)
# cd $dbs
# mkdir -p COSMIC
# cd COSMIC
# # manual download of CancerMutationCensus_AllData_Tsv_v99_GRCh38.tar, Cosmic_GenomeScreensMutant_Vcf_v99_GRCh38.tar, Cosmic_CompleteTargetedScreensMutant_Vcf_v99_GRCh38.tar and Cosmic_NonCodingVariants_Vcf_v99_GRCh38.tar from https://apps.ingenuity.com/ingsso/login
# ls *.tar | xargs -l1 tar -xf 
# gunzip -c CancerMutationCensus_AllData_v99_GRCh38.tsv.gz | php $src/Install/db_converter_cosmic.php -in_genome_vcf Cosmic_GenomeScreensMutant_v99_GRCh38.vcf.gz -in_non_coding_vcf Cosmic_NonCodingVariants_v99_GRCh38.vcf.gz -in_target_screens_vcf Cosmic_CompleteTargetedScreensMutant_v99_GRCh38.vcf.gz -out cmc_export_v99.vcf.gz

# # when using the containerized megSAP version:
# cd <path-to-host-data-folder>/dbs/
# mkdir -p COSMIC
# cd COSMIC
# # manual download of CancerMutationCensus_AllData_Tsv_v99_GRCh38.tar, Cosmic_GenomeScreensMutant_Vcf_v99_GRCh38.tar, Cosmic_CompleteTargetedScreensMutant_Vcf_v99_GRCh38.tar and Cosmic_NonCodingVariants_Vcf_v99_GRCh38.tar from https://apps.ingenuity.com/ingsso/login
# ls *.tar | xargs -l1 tar -xf 
#apptainer exec -B <path-to-host-data-folder>:/megSAP/data/data_folder/ --pwd /megSAP/data/data_folder/dbs/COSMIC megSAP_[version].sif sh -c "gunzip -c CancerMutationCensus_AllData_v99_GRCh38.tsv.gz | php /megSAP/src/Install/db_converter_cosmic.php -in_genome_vcf Cosmic_GenomeScreensMutant_v99_GRCh38.vcf.gz -in_non_coding_vcf Cosmic_NonCodingVariants_v99_GRCh38.vcf.gz -in_target_screens_vcf Cosmic_CompleteTargetedScreensMutant_v99_GRCh38.vcf.gz -out cmc_export_v99.vcf.gz"