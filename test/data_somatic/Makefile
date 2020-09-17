all: clear link case1 case2 case3

clear:
	rm -rf Sample_GS110168_03 Sample_GS110169_03 Sample_GS150344_01
	rm -rf out_case1 out_case2 out_case3

link:
	mkdir -p Sample_GS110168_03 Sample_GS110169_03 Sample_GS150344_01 Sample_RX000001_01
	mkdir -p out_case1 out_case2 out_case3
	# input files for case 1
	php ../link_test_data.php GS150344_01_GCTCGGTAextracted_L001_R1_001.fastq.gz Sample_GS150344_01/GS150344_01_GCTCGGTAextracted_L001_R1_001.fastq.gz
	php ../link_test_data.php GS150344_01_GCTCGGTAextracted_L001_R2_001.fastq.gz Sample_GS150344_01/GS150344_01_GCTCGGTAextracted_L001_R2_001.fastq.gz
	php ../link_test_data.php GS150344_01_GCTCGGTAextracted_L001_R1_001.fastq.gz Sample_GS150344_01/GS150344_01_GCTCGGTAextracted_L001_R1_001.fastq.gz
	php ../link_test_data.php GS150344_01_GCTCGGTAextracted_L001_R2_001.fastq.gz Sample_GS150344_01/GS150344_01_GCTCGGTAextracted_L001_R2_001.fastq.gz
	# input files for case 2
	php ../link_test_data.php GS110168_chr21_L007_R1_001.fastq.gz Sample_GS110168_03/GS110168_03_chr21_L007_R1_001.fastq.gz
	php ../link_test_data.php GS110168_chr21_L007_R2_001.fastq.gz Sample_GS110168_03/GS110168_03_chr21_L007_R2_001.fastq.gz
	php ../link_test_data.php GS110169_chr21_L008_R1_001.fastq.gz Sample_GS110169_03/GS110169_03_chr21_L008_R1_001.fastq.gz
	php ../link_test_data.php GS110169_chr21_L008_R2_001.fastq.gz Sample_GS110169_03/GS110169_03_chr21_L008_R2_001.fastq.gz
	# input RNA data for an_rna_step of case 2
	cp ../data/somatic_pipeline_rna_counts.tsv Sample_RX000001_01/RX000001_01_counts.tsv
	cp ../data/somatic_pipeline_tumor_rna.bam Sample_RX000001_01/RX000001_01.bam
	cp ../data/somatic_pipeline_tumor_rna.bam.bai Sample_RX000001_01/RX000001_01.bam.bai

# tumor-only
case1: case1_analyze case1_somatic case1_diff

case1_analyze:
	php ../../src/Pipelines/analyze.php \
		-folder Sample_GS150344_01 \
		-name GS150344_01 \
		-steps ma -somatic \
		-threads 1 \
		-system system_single.ini \
		--log Sample_GS150344_01/analyze.log

case1_somatic:
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
		-folder Sample_GS110168_03 \
		-name GS110168_03 \
		-steps ma -somatic \
		-system system_pair.ini \
		-threads 1 \
		--log Sample_GS110168_03/analyze.log
	php ../../src/Pipelines/analyze.php \
		-folder Sample_GS110169_03 \
		-name GS110169_03 \
		-steps ma,vc,cn -somatic \
		-system system_pair.ini \
		-threads 1 \
		--log Sample_GS110169_03/analyze.log
case2_somatic:
	php ../../src/Pipelines/somatic_dna.php \
		-t_bam Sample_GS110168_03/GS110168_03.bam -n_bam Sample_GS110169_03/GS110169_03.bam \
		-out_folder out_case2 \
		-prefix GS110168_03-GS110169_03 \
		-steps vc,an,msi,an_rna,ci\
		-t_rna_bam Sample_RX000001_01/RX000001_01.bam \
		-rna_ref_tissue gallbladder \
		-system system_pair.ini -n_system system_pair.ini \
		-cancer_type CANCER \
		-include_germline \
		--log out_case2/GS110168_03-GS110169_03.log
case2_diff:
	# variants
	php ../compare_variants.php -in1 expected_case2/GS110168_03-GS110169_03.GSvar -in2 out_case2/GS110168_03-GS110169_03.GSvar
	php ../compare_bed.php expected_case2/GS110168_03-GS110169_03_stat_lowcov.bed out_case2/GS110168_03-GS110169_03_stat_lowcov.bed
	test -s out_case2/GS110168_03-GS110169_03_cgi_drug_prescription.tsv
	test -s out_case2/GS110168_03-GS110169_03_cgi_drug_prescription_bioactivities.tsv
	test -s out_case2/GS110168_03-GS110169_03_cgi_mutation_analysis.tsv
	test -s Sample_GS110169_03/GS110169_03_cgi_mutation_analysis.tsv
	# qcML files
	php ../compare_qcml_xml.php expected_case2/GS110168_03-GS110169_03_stats_som.qcML out_case2/GS110168_03-GS110169_03_stats_som.qcML
	# MSI file
	numdiff expected_case2/GS110168_03-GS110169_03_msi.tsv out_case2/GS110168_03-GS110169_03_msi.tsv --absolute-tolerance 0.01
	# bafs
	numdiff expected_case2/GS110168_03-GS110169_03_bafs.igv out_case2/GS110168_03-GS110169_03_bafs.igv --absolute-tolerance 0.02

#Test MANTA SV calls with positive fusion control
case3: case3_somatic case3_diff
case3_somatic:
	php ../../src/Pipelines/somatic_dna.php \
		-t_bam ../data/manta_sv_tumor.bam -n_bam ../data/manta_sv_normal.bam \
		-out_folder out_case3 \
		-prefix manta_sv_tumor-manta_sv_normal \
		-steps vc\
		-system system_pair_sv.ini -n_system system_pair_sv.ini \
		-skip_correlation \
		--log out_case3/tumor-normal.log
case3_diff:
	diff expected_case3/manta_sv_tumor-manta_sv_normal_manta_var_structural.bedpe out_case3/manta_sv_tumor-manta_sv_normal_manta_var_structural.bedpe --ignore-matching-lines="##cmdline" --ignore-matching-lines="##fileDate="
