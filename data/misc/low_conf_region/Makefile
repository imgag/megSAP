help:
	@cat Makefile

########################## error regions based on trios ##########################

trios:
	php trios.php > trios.tsv

trio_errors:
	php trio_errors.php trios.tsv ssHAEv6 4 > trio_errors.tsv
	php trio_errors.php trios.tsv ssHAEv7 3 >> trio_errors.tsv
	php trio_errors.php trios.tsv TruSeqPCRfree 8 >> trio_errors.tsv	
	php gsvar2bed.php -in trio_errors_ssHAEv6/*.GSvar -out trio_errors_ssHAEv6.bed -errors 2 > trio_errors_ssHAEv6.log
	php gsvar2bed.php -in trio_errors_ssHAEv7/*.GSvar -out trio_errors_ssHAEv7.bed -errors 2 > trio_errors_ssHAEv7.log
	php gsvar2bed.php -in trio_errors_TruSeqPCRfree/*.GSvar -out trio_errors_TruSeqPCRfree.bed -errors 2 > trio_errors_TruSeqPCRfree.log
	BedAdd -in trio_errors_ssHAEv6.bed trio_errors_ssHAEv7.bed trio_errors_TruSeqPCRfree.bed | BedMerge -out trio_errors.bed
	BedExtend -n 5 -in trio_errors.bed | BedMerge | BedShrink -n 5 >  trio_errors_moc5.bed
	BedSubtract -in /mnt/storage2/GRCh38/share/data/enrichment/WGS_hg38.bed -in2 trio_errors_moc5.bed -out trio_errors_moc5_high_conf.bed

########################## error regions based on twins ##########################

twins:
	php twins.php > twins.tsv

twin_vars_discordant:
	mkdir -p twins
	php vars_discordant.php twins.tsv

twin_errors:
	php twin_errors.php twin_errors.bed twins/*_disc.vcf
	BedExtend -n 5 -in twin_errors.bed | BedMerge | BedShrink -n 5 >  twin_errors_moc5.bed
	BedSubtract -in /mnt/storage2/GRCh38/share/data/enrichment/WGS_hg38.bed -in2 twin_errors_moc5.bed -out twin_errors_moc5_high_conf.bed

########################## gnomAD AC=0 ##########################

gnomad:
	php gnomad.php /mnt/storage2/GRCh38/share/data/dbs/gnomAD/gnomAD_genome_v3.1.1_GRCh38.vcf.gz AC0 4 gnomAD_filter.bed
	BedExtend -n 5 -in gnomAD_filter.bed | BedMerge | BedShrink -n 6 | BedExtend -n 1  > gnomAD_filter_moc5.bed
	BedSubtract -in /mnt/storage2/GRCh38/share/data/enrichment/WGS_hg38.bed -in2 gnomAD_filter_moc5.bed -out gnomAD_filter_moc5_high_conf.bed

########################## combined above low-confidence regions ##########################

info:
	BedInfo -filename -in twin_errors_moc5.bed | grep Bases
	BedInfo -filename -in trio_errors_moc5.bed | grep Bases
	BedInfo -filename -in gnomAD_filter_moc5.bed | grep Bases

combined:
	BedAdd -in trio_errors_moc5.bed twin_errors_moc5.bed gnomAD_filter_moc5.bed | BedMerge -out combined.bed
	BedExtend -n 5 -in combined.bed | BedMerge | BedShrink -n 5 >  combined_moc5.bed
	BedSubtract -in /mnt/storage2/GRCh38/share/data/enrichment/WGS_hg38.bed -in2 combined_moc5.bed -out combined_moc5_high_conf.bed
	BedExtend -n 10 -in combined.bed | BedMerge | BedShrink -n 10 >  combined_moc10.bed
	BedSubtract -in /mnt/storage2/GRCh38/share/data/enrichment/WGS_hg38.bed -in2 combined_moc10.bed -out combined_moc10_high_conf.bed
	BedExtend -n 20 -in combined.bed | BedMerge | BedShrink -n 20 >  combined_moc20.bed
	BedSubtract -in /mnt/storage2/GRCh38/share/data/enrichment/WGS_hg38.bed -in2 combined_moc20.bed -out combined_moc20_high_conf.bed
	BedExtend -n 50 -in combined.bed | BedMerge | BedShrink -n 50 >  combined_moc50.bed
	BedSubtract -in /mnt/storage2/GRCh38/share/data/enrichment/WGS_hg38.bed -in2 combined_moc50.bed -out combined_moc50_high_conf.bed
	BedExtend -n 1 -in combined_moc50.bed | BedMerge > combined_moc50_ext1.bed
	BedSubtract -in /mnt/storage2/GRCh38/share/data/enrichment/WGS_hg38.bed -in2 combined_moc50_ext1.bed -out combined_moc50_ext1_high_conf.bed
	BedExtend -n 2 -in combined_moc50.bed | BedMerge > combined_moc50_ext2.bed
	BedSubtract -in /mnt/storage2/GRCh38/share/data/enrichment/WGS_hg38.bed -in2 combined_moc50_ext2.bed -out combined_moc50_ext2_high_conf.bed
	BedExtend -n 5 -in combined_moc50.bed | BedMerge > combined_moc50_ext5.bed
	BedSubtract -in /mnt/storage2/GRCh38/share/data/enrichment/WGS_hg38.bed -in2 combined_moc50_ext5.bed -out combined_moc50_ext5_high_conf.bed
	BedExtend -n 10 -in combined_moc50.bed | BedMerge > combined_moc50_ext10.bed
	BedSubtract -in /mnt/storage2/GRCh38/share/data/enrichment/WGS_hg38.bed -in2 combined_moc50_ext10.bed -out combined_moc50_ext10_high_conf.bed

########################## benchmarking ##########################

gene_exons:
	NGSDExportGenes | grep "protein-coding gene" | cut -f1 | GenesToBed -source ensembl -mode exon | BedExtend -n 5 | BedMerge > gene_exons.bed

patho:
	zcat /mnt/storage2/GRCh38/share/data/dbs/NGSD/NGSD_germline.vcf.gz | egrep "^#|CLAS=4|CLAS=5" | grep -v chrY | grep -v chrMT > NGSD_patho.vcf
	zcat /mnt/storage2/GRCh38/share/data/dbs/ClinVar/clinvar_20211212_converted_GRCh38.vcf.gz | egrep "^#|pathogenic" | grep -v conflicting_interpretations_of_pathogenicity | grep -v chrY | grep -v chrMT  > ClinVar_patho.vcf
	php ngsd_variants.php causal > NGSD_causal.vcf
	php ngsd_variants.php artefact > NGSD_artefact.vcf

benchmark:
	php benchmark.php twins/*_disc.vcf NGSD_artefact.vcf NGSD_causal.vcf NGSD_patho.vcf ClinVar_patho.vcf /mnt/storage2/GRCh38/share/data/dbs/GIAB/NA12878/high_conf_regions.bed /mnt/storage2/GRCh38/share/data/dbs/GIAB/NA12877/high_conf_regions.bed gnomAD_filter_moc5_high_conf.bed trio_errors_moc5_high_conf.bed twin_errors_moc5_high_conf.bed combined_moc5_high_conf.bed combined_moc10_high_conf.bed combined_moc20_high_conf.bed combined_moc50_high_conf.bed combined_moc50_ext1_high_conf.bed combined_moc50_ext2_high_conf.bed combined_moc50_ext5_high_conf.bed combined_moc50_ext10_high_conf.bed> benchmark.tsv

########################## deploy ##########################

deploy:
	cp combined_moc50_ext1.bed ../low_conf_regions.bed
	
	
