all: clear link case1 case2 case4

clear:
	rm -rf Sample_GS110168 Sample_GS110169 Sample_GS150344_01 Sample_Tumor2 Sample_Normal2
	rm -rf out_case1 out_case2 out_case3 out_case4

link:
	mkdir -p Sample_GS110169 Sample_GS110168 Sample_GS150344_01 Sample_GS150344_01 Sample_Tumor2 Sample_Normal2
	mkdir -p out_case1 out_case2 out_case3 out_case4
	# input files for case 1
	php ../link_test_data.php GS150344_01_GCTCGGTAextracted_L001_R1_001.fastq.gz Sample_GS150344_01/GS150344_01_GCTCGGTAextracted_L001_R1_001.fastq.gz
	php ../link_test_data.php GS150344_01_GCTCGGTAextracted_L001_R2_001.fastq.gz Sample_GS150344_01/GS150344_01_GCTCGGTAextracted_L001_R2_001.fastq.gz
	php ../link_test_data.php GS150344_01_GCTCGGTAextracted_L001_R1_001.fastq.gz Sample_GS150344_01/GS150344_01_GCTCGGTAextracted_L001_R1_001.fastq.gz
	php ../link_test_data.php GS150344_01_GCTCGGTAextracted_L001_R2_001.fastq.gz Sample_GS150344_01/GS150344_01_GCTCGGTAextracted_L001_R2_001.fastq.gz
	# input files for case 2, 3
	php ../link_test_data.php GS110168_chr21_L007_R1_001.fastq.gz Sample_GS110168/GS110168_chr21_L007_R1_001.fastq.gz
	php ../link_test_data.php GS110168_chr21_L007_R2_001.fastq.gz Sample_GS110168/GS110168_chr21_L007_R2_001.fastq.gz
	php ../link_test_data.php GS110169_chr21_L008_R1_001.fastq.gz Sample_GS110169/GS110169_chr21_L008_R1_001.fastq.gz
	php ../link_test_data.php GS110169_chr21_L008_R2_001.fastq.gz Sample_GS110169/GS110169_chr21_L008_R2_001.fastq.gz
	# input files for case 4
	php ../link_test_data.php Tumor2_XX_L00X_R1_001.fastq.gz Sample_Tumor2/Tumor2_XX_L00X_R1_001.fastq.gz
	php ../link_test_data.php Tumor2_XX_L00X_R2_001.fastq.gz Sample_Tumor2/Tumor2_XX_L00X_R2_001.fastq.gz
	php ../link_test_data.php Normal2_XX_L00X_R1_001.fastq.gz Sample_Normal2/Normal2_XX_L00X_R1_001.fastq.gz
	php ../link_test_data.php Normal2_XX_L00X_R2_001.fastq.gz Sample_Normal2/Normal2_XX_L00X_R2_001.fastq.gz

# tumor-only
case1: case1_run case1_diff

case1_run:
	php ../../src/Pipelines/analyze.php \
		-folder Sample_GS150344_01 \
		-name GS150344_01 \
		-steps ma -somatic \
		-system system_single.ini \
		--log Sample_GS150344_01/analyze.log
	php ../../src/Pipelines/somatic_dna.php \
		-t_bam Sample_GS150344_01/GS150344_01.bam \
		-out_folder out_case1 \
		-prefix GS150344_01 \
		-steps vc,cn,an \
		-system system_single.ini \
		--log out_case1/GS150344.log

case1_diff:
	php ../compare_bed.php expected_case1/GS150344_01_stat_lowcov.bed out_case1/GS150344_01_stat_lowcov.bed
	# variants
	php ../compare_variants.php -in1 expected_case1/GS150344_01.GSvar -in2 out_case1/GS150344_01.GSvar
	test -s out_case1/GS150344_01_cnvs.tsv

# tumor-normal (strelka)
case2: case2_analyze case2_somatic case2_diff

case2_analyze:
	php ../../src/Pipelines/analyze.php \
		-folder Sample_GS110168 \
		-name GS110168 \
		-steps ma -somatic \
		-system system_pair.ini \
		--log Sample_GS110168/analyze.log
	php ../../src/Pipelines/analyze.php \
		-folder Sample_GS110169 \
		-name GS110169 \
		-steps ma -somatic \
		-system system_pair.ini \
		--log Sample_GS110169/analyze.log

case2_somatic:
	php ../../src/Pipelines/somatic_dna.php \
		-t_bam Sample_GS110168/GS110168.bam -n_bam Sample_GS110169/GS110169.bam \
		-out_folder out_case2 \
		-prefix GS110168-GS110169 \
		-steps vc,cn,an,ci,msi \
		-cancer_type CANCER \
		-system system_pair.ini \
		-n_system system_pair.ini \
		--log out_case2/GS110168-GS110169.log

case2_diff:
	# variants
	php ../compare_variants.php -in1 expected_case2/GS110168-GS110169.GSvar -in2 out_case2/GS110168-GS110169.GSvar
	php ../compare_variants.php -in1 expected_case2/GS110168-GS110169_var_qci.vcf.gz -in2 out_case2/GS110168-GS110169_var_qci.vcf.gz
	php ../compare_bed.php expected_case2/GS110168-GS110169_stat_lowcov.bed out_case2/GS110168-GS110169_stat_lowcov.bed
	
	test -s out_case2/GS110168-GS110169_cgi_drug_prescription.tsv
	test -s out_case2/GS110168-GS110169_cgi_drug_prescription_bioactivities.tsv
	test -s out_case2/GS110168-GS110169_cgi_mutation_analysis.tsv
	test -s out_case2/GS110168-GS110169_cnvs.tsv

	# qcML files
	php ../compare_qcml_xml.php expected_case2/GS110168-GS110169_stats_som.qcML out_case2/GS110168-GS110169_stats_som.qcML

	# bafs
	diff expected_case2/GS110168-GS110169_bafs.igv out_case2/GS110168-GS110169_bafs.igv
	
	# MSI file
	diff expected_case2/GS110168-GS110169_msi.tsv out_case2/GS110168-GS110169_msi.tsv

# tumor-normal (strelka)
case4: case4_analyze case4_somatic case4_diff

case4_analyze:
	php ../../src/Pipelines/analyze.php \
		-folder Sample_Tumor2 \
		-name Tumor2 \
		-steps ma -somatic \
		-system system_pair.ini \
		--log Sample_Tumor2/analyze.log
	php ../../src/Pipelines/analyze.php \
		-folder Sample_Normal2 \
		-name Normal2 \
		-steps ma,vc -somatic \
		-system system_pair.ini \
		--log Sample_Normal2/analyze.log

case4_somatic:
	php ../../src/Pipelines/somatic_dna.php \
		-t_bam Sample_Tumor2/Tumor2.bam -n_bam Sample_Normal2/Normal2.bam \
		-out_folder out_case4 \
		-prefix Tumor2-Normal2 \
		-steps vc,cn,an,ci \
		-system system_pair.ini -n_system system_pair.ini \
		-include_germline \
		-cancer_type CANCER \
		--log out_case4/Tumor2-Normal2.log

case4_diff:
	# variants
	php ../compare_variants.php -in1 expected_case4/Tumor2-Normal2.GSvar -in2 out_case4/Tumor2-Normal2.GSvar
	php ../compare_variants.php -in1 expected_case4/Tumor2-Normal2_var_qci.vcf.gz -in2 out_case4/Tumor2-Normal2_var_qci.vcf.gz
	test -s out_case4/Tumor2-Normal2_cgi_drug_prescription.tsv
	test -s out_case4/Tumor2-Normal2_cgi_drug_prescription_bioactivities.tsv
	test -s out_case4/Tumor2-Normal2_cgi_mutation_analysis.tsv
	test -s out_case4/Tumor2-Normal2_cnvs.tsv

	# qcML files
	php ../compare_qcml_xml.php expected_case4/Tumor2-Normal2_stats_som.qcML out_case4/Tumor2-Normal2_stats_som.qcML
	
	# bafs
	diff expected_case4/Tumor2-Normal2_bafs.igv out_case4/Tumor2-Normal2_bafs.igv
