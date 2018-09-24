all: clear link case1 case2 case4

clear:
	rm -rf Sample_GS110168 Sample_GS110169 Sample_GS150344_01  Sample_DX183978_01 Sample_DX183927_01
	rm -rf out_case1 out_case2 out_case3

link:
	mkdir -p Sample_GS110168 Sample_GS110169 Sample_GS150344_01 Sample_DX183978_01 Sample_DX183927_01
	mkdir -p out_case1 out_case2 out_case3
	# input files for case 1
	php ../link_test_data.php GS150344_01_GCTCGGTAextracted_L001_R1_001.fastq.gz Sample_GS150344_01/GS150344_01_GCTCGGTAextracted_L001_R1_001.fastq.gz
	php ../link_test_data.php GS150344_01_GCTCGGTAextracted_L001_R2_001.fastq.gz Sample_GS150344_01/GS150344_01_GCTCGGTAextracted_L001_R2_001.fastq.gz
	php ../link_test_data.php GS150344_01_GCTCGGTAextracted_L001_R1_001.fastq.gz Sample_GS150344_01/GS150344_01_GCTCGGTAextracted_L001_R1_001.fastq.gz
	php ../link_test_data.php GS150344_01_GCTCGGTAextracted_L001_R2_001.fastq.gz Sample_GS150344_01/GS150344_01_GCTCGGTAextracted_L001_R2_001.fastq.gz
	# input files for case 2
	php ../link_test_data.php GS110168_chr21_L007_R1_001.fastq.gz Sample_GS110168/GS110168_chr21_L007_R1_001.fastq.gz
	php ../link_test_data.php GS110168_chr21_L007_R2_001.fastq.gz Sample_GS110168/GS110168_chr21_L007_R2_001.fastq.gz
	php ../link_test_data.php GS110169_chr21_L008_R1_001.fastq.gz Sample_GS110169/GS110169_chr21_L008_R1_001.fastq.gz
	php ../link_test_data.php GS110169_chr21_L008_R2_001.fastq.gz Sample_GS110169/GS110169_chr21_L008_R2_001.fastq.gz
	#input files for case 3
	php ../link_test_data.php DX183978_01_chr1_chr2.bam Sample_DX183978_01/DX183978_01.bam
	php ../link_test_data.php DX183978_01_chr1_chr2.bam.bai Sample_DX183978_01/DX183978_01.bam.bai
	php ../link_test_data.php DX183927_01_chr1_chr2.bam Sample_DX183927_01/DX183927_01.bam
	php ../link_test_data.php DX183927_01_chr1_chr2.bam.bai Sample_DX183927_01/DX183927_01.bam.bai

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
		-steps vc,an,cn \
		-system system_single.ini \
		--log out_case1/GS150344.log
case1_diff:
	php ../compare_bed.php expected_case1/GS150344_01_stat_lowcov.bed out_case1/GS150344_01_stat_lowcov.bed
	php ../compare_variants.php -in1 expected_case1/GS150344_01.GSvar -in2 out_case1/GS150344_01.GSvar
	test -s out_case1/GS150344_01_cnvs.tsv

# tumor-normal (strelka), testing without ClinCNV (GC content too low in input files reduced to one chr)
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
		-steps vc,an,msi,ci\
		-system system_pair.ini -n_system system_pair.ini \
		-cancer_type CANCER \
		--log out_case2/GS110168-GS110169.log
case2_diff:
	# variants
	php ../compare_variants.php -in1 expected_case2/GS110168-GS110169.GSvar -in2 out_case2/GS110168-GS110169.GSvar
	php ../compare_variants.php -in1 expected_case2/GS110168-GS110169_var_qci.vcf.gz -in2 out_case2/GS110168-GS110169_var_qci.vcf.gz
	php ../compare_bed.php expected_case2/GS110168-GS110169_stat_lowcov.bed out_case2/GS110168-GS110169_stat_lowcov.bed
	test -s out_case2/GS110168-GS110169_cgi_drug_prescription.tsv
	test -s out_case2/GS110168-GS110169_cgi_drug_prescription_bioactivities.tsv
	test -s out_case2/GS110168-GS110169_cgi_mutation_analysis.tsv
	# qcML files
	php ../compare_qcml_xml.php expected_case2/GS110168-GS110169_stats_som.qcML out_case2/GS110168-GS110169_stats_som.qcML
	# bafs
	diff expected_case2/GS110168-GS110169_bafs.igv out_case2/GS110168-GS110169_bafs.igv
	# MSI file
	diff expected_case2/GS110168-GS110169_msi.tsv out_case2/GS110168-GS110169_msi.tsv

# tumor-normal, intended for testing CNV calls, bam file already exists
case3: case3_somatic case3_diff
case3_somatic:
	php ../../src/Pipelines/somatic_dna.php \
		-t_bam Sample_DX183978_01/DX183978_01.bam -n_bam Sample_DX183927_01/DX183927_01.bam \
		-out_folder out_case3 \
		-prefix DX183978_01-DX183927_01 \
		-steps vc,an,cn,ci\
		-system somatic_cn_test_pair.ini -n_system somatic_cn_test_pair.ini \
		-cancer_type CANCER \
		--log out_case3/DX183978_01-DX183927_01.log
case3_diff:
	diff expected_case3/DX183978_01-DX183927_01_clincnv.tsv out_case3/DX183978_01-DX183927_01_clincnv.tsv
	diff expected_case3/DX183978_01-DX183927_01.GSvar out_case3/DX183978_01-DX183927_01.GSvar
