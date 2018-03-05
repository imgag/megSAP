all: clear link case1 case2 case3 case4

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
	php ../../src/Pipelines/somatic_dna.php \
	    -p_folder . \
	    -t_id GS150344_01 -n_id na \
	    -o_folder out_case1 \
	    -steps ma,vc,an \
	    -t_sys system_single.ini \
	    --log out_case1/GS150344.log

case1_diff:
	php ../compare_bed.php expected_case1/GS150344_01_stat_lowcov.bed out_case1/GS150344_01_stat_lowcov.bed
	# variants
	php ../compare_variants.php expected_case1/GS150344_01.GSvar out_case1/GS150344_01.GSvar
	test -s out_case1/GS150344_01_cnvs.tsv

# tumor-normal (strelka)
case2: case2_run case2_diff

case2_run:
	php ../../src/Pipelines/somatic_dna.php \
	    -p_folder . \
	    -t_id GS110168 -n_id GS110169 \
	    -o_folder out_case2 \
	    -steps ma,vc,an,ci \
	    -reduce_variants_filter \
	    -t_sys system_pair.ini -n_sys system_pair.ini \
	    --log out_case2/GS110168-GS110169.log

case2_diff:
	# variants
	php ../compare_variants.php expected_case2/GS110168-GS110169.GSvar out_case2/GS110168-GS110169.GSvar
	php ../compare_variants.php expected_case2/GS110168-GS110169_var_qci.vcf.gz out_case2/GS110168-GS110169_var_qci.vcf.gz
	php ../compare_bed.php expected_case2/GS110168-GS110169_stat_lowcov.bed out_case2/GS110168-GS110169_stat_lowcov.bed
	
	test -s out_case2/GS110168-GS110169_cgi_drug_prescription.tsv
	test -s out_case2/GS110168-GS110169_cgi_drug_prescription_bioactivities.tsv
	test -s out_case2/GS110168-GS110169_cgi_mutation_analysis.tsv
	test -s out_case2/GS110168-GS110169_cnvs.tsv

	# qcML files
	php ../compare_qcml_xml.php expected_case2/GS110168-GS110169_stats_som.qcML out_case2/GS110168-GS110169_stats_som.qcML

	# bafs
	diff expected_case2/GS110168-GS110169_bafs.igv out_case2/GS110168-GS110169_bafs.igv

# tumor-normal (freebayes)
# depends on case2_run
case3: case3_run case3_diff

case3_run:
	php ../../src/Pipelines/somatic_dna.php \
	    -p_folder . \
	    -t_id GS110168 -n_id GS110169 \
	    -o_folder out_case3 \
	    -steps vc,an,ci \
	    -t_sys system_pair.ini -n_sys system_pair.ini \
	    -freebayes \
	    --log out_case3/GS110168-GS110169.log

case3_diff:
	# variants
	php ../compare_variants.php expected_case3/GS110168-GS110169.GSvar out_case3/GS110168-GS110169.GSvar
	php ../compare_variants.php expected_case3/GS110168-GS110169_var_qci.vcf.gz out_case3/GS110168-GS110169_var_qci.vcf.gz
	test -s out_case3/GS110168-GS110169_cgi_drug_prescription.tsv
	test -s out_case3/GS110168-GS110169_cgi_drug_prescription_bioactivities.tsv
	test -s out_case3/GS110168-GS110169_cgi_mutation_analysis.tsv
	test -s out_case3/GS110168-GS110169_cnvs.tsv
	
	# qcML files
	php ../compare_qcml_xml.php expected_case3/GS110168-GS110169_stats_som.qcML out_case3/GS110168-GS110169_stats_som.qcML

	# bafs
	diff expected_case3/GS110168-GS110169_bafs.igv out_case3/GS110168-GS110169_bafs.igv

# tumor-normal (strelka)
case4: case4_run case4_diff

case4_run:
	php ../../src/Pipelines/somatic_dna.php \
	    -p_folder . \
	    -t_id Tumor2 -n_id Normal2 \
	    -o_folder out_case4 \
	    -steps ma,vc,an,ci \
	    -t_sys system_pair.ini -n_sys system_pair.ini \
	    --log out_case4/Tumor2-Normal2.log
	
case4_diff:
	# variants
	php ../compare_variants.php expected_case4/Tumor2-Normal2.GSvar out_case4/Tumor2-Normal2.GSvar
	php ../compare_variants.php expected_case4/Tumor2-Normal2_var_qci.vcf.gz out_case4/Tumor2-Normal2_var_qci.vcf.gz
	test -s out_case4/Tumor2-Normal2_cgi_drug_prescription.tsv
	test -s out_case4/Tumor2-Normal2_cgi_drug_prescription_bioactivities.tsv
	test -s out_case4/Tumor2-Normal2_cgi_mutation_analysis.tsv
	test -s out_case4/Tumor2-Normal2_cnvs.tsv

	# qcML files
	php ../compare_qcml_xml.php expected_case4/Tumor2-Normal2_stats_som.qcML out_case4/Tumor2-Normal2_stats_som.qcML
	
	# bafs
	diff expected_case4/Tumor2-Normal2_bafs.igv out_case4/Tumor2-Normal2_bafs.igv