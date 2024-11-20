#!/bin/bash
set -e
set -o pipefail
set -o verbose

root=`pwd`
ngsbits=$root/tools/ngs-bits/bin/
genome=$root/genomes/GRCh38.fa
giab=$root/dbs/GIAB/

#create GIAB folder if missing
mkdir -p $giab/

#download and normalize HG001/NA12878
mkdir $giab/NA12878/
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -O $giab/NA12878/high_conf_variants.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.bed -O $giab/NA12878/high_conf_regions.bed
zcat $giab/NA12878/high_conf_variants.vcf.gz | $ngsbits/VcfBreakMulti | $ngsbits/VcfFilter -remove_invalid | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | bgzip > $giab/NA12878/high_conf_variants_normalized.vcf.gz
tabix $giab/NA12878/high_conf_variants_normalized.vcf.gz

#download and normalize HG002/NA24385
mkdir $giab/NA24385/
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -O $giab/NA24385/high_conf_variants.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed -O $giab/NA24385/high_conf_regions.bed
zcat $giab/NA24385/high_conf_variants.vcf.gz | $ngsbits/VcfBreakMulti | $ngsbits/VcfFilter -remove_invalid | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | bgzip > $giab/NA24385/high_conf_variants_normalized.vcf.gz
tabix $giab/NA24385/high_conf_variants_normalized.vcf.gz

#download and normalize HG002/NA24385 CMRG
mkdir $giab/NA24385_CMRG/
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/SmallVariant/HG002_GRCh38_CMRG_smallvar_v1.00.vcf.gz -O $giab/NA24385_CMRG/high_conf_variants.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/SmallVariant/HG002_GRCh38_CMRG_smallvar_v1.00.bed -O $giab/NA24385_CMRG/high_conf_regions.bed
zcat $giab/NA24385_CMRG/high_conf_variants.vcf.gz | $ngsbits/VcfBreakMulti | $ngsbits/VcfFilter -remove_invalid | $ngsbits/VcfLeftNormalize -stream -ref $genome | $ngsbits/VcfStreamSort | bgzip > $giab/NA24385_CMRG/high_conf_variants_normalized.vcf.gz
tabix $giab/NA24385_CMRG/high_conf_variants_normalized.vcf.gz
