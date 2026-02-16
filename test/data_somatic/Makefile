SHELL=/bin/bash
all: clear link tumor_only tumor_normal tumor_normal_with_fusion

clear:
	rm -rf Sample_GS110168_03 Sample_GS110169_03 Sample_DX000001_01 Sample_RX000001_01
	rm -rf out_case1 out_case2 out_case3

link:
	mkdir -p Sample_GS110168_03 Sample_GS110169_03 Sample_DX000001_01 Sample_RX000001_01
	mkdir -p out_case1 out_case2 out_case3
	# input files for case 1
	php ../link_test_data.php som_pipeline_test_case1_L001_R1_001.fastq.gz Sample_DX000001_01/DX000001_01_L001_R1_001.fastq.gz
	php ../link_test_data.php som_pipeline_test_case1_L001_R2_001.fastq.gz Sample_DX000001_01/DX000001_01_L001_R2_001.fastq.gz
	# input files for case 2
	php ../link_test_data.php GS110168_chr21_L007_R1_001.fastq.gz Sample_GS110168_03/GS110168_03_chr21_L007_R1_001.fastq.gz
	php ../link_test_data.php GS110168_chr21_L007_R2_001.fastq.gz Sample_GS110168_03/GS110168_03_chr21_L007_R2_001.fastq.gz
	php ../link_test_data.php GS110169_chr21_L008_R1_001.fastq.gz Sample_GS110169_03/GS110169_03_chr21_L008_R1_001.fastq.gz
	php ../link_test_data.php GS110169_chr21_L008_R2_001.fastq.gz Sample_GS110169_03/GS110169_03_chr21_L008_R2_001.fastq.gz
	# input RNA data for an_rna_step of case 2
	cp ../data/somatic_pipeline_rna_counts.tsv Sample_RX000001_01/RX000001_01_counts.tsv
	cp ../data/somatic_pipeline_tumor_rna.bam Sample_RX000001_01/RX000001_01.bam
	cp ../data/somatic_pipeline_tumor_rna.bam.bai Sample_RX000001_01/RX000001_01.bam.bai

# tumor-only (varscan) - without CNVs
tumor_only: to_analyze to_somatic to_diff

to_analyze:
	php ../../src/Pipelines/analyze.php \
		-folder Sample_DX000001_01 \
		-name DX000001_01 \
		-steps ma -somatic \
		-threads 1 \
		-system system_single.ini \
		--log Sample_DX000001_01/analyze.log

to_somatic:
	php ../../src/Pipelines/somatic_tumor_only.php \
		-t_bam Sample_DX000001_01/DX000001_01.cram \
		-out_folder out_case1 \
		-prefix DX000001_01 \
		-steps vc,vi,an \
		-system system_single.ini \
		-min_depth_t 100 \
		--log out_case1/DX000001_01.log

to_diff:
	# low coverage bed
	php ../compare_bed.php expected_case1/DX000001_01_stat_lowcov.bed out_case1/DX000001_01_stat_lowcov.bed
	# small variants
	php ../compare_variants.php -in1 expected_case1/DX000001_01.GSvar -in2 out_case1/DX000001_01.GSvar -add tumor_af,tumor_dp
	# virus detection
	diff expected_case1/DX000001_01_viral.tsv Sample_DX000001_01/DX000001_01_viral.tsv
	# qcML files
	php ../compare_qcml_xml.php expected_case1/DX000001_01_stats_other.qcML out_case1/DX000001_01_stats_other.qcML

# tumor-normal (strelka) - without CNVs (GC content too low in input files reduced to one chr)
tumor_normal: tn_analyze tn_somatic tn_diff

tn_analyze:
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
		-steps ma,vc -somatic \
		-system system_pair.ini \
		-threads 1 \
		--log Sample_GS110169_03/analyze.log

tn_somatic:
	rm -rf out_case2/GS110168_03-GS110169_03.log
	php ../../src/Pipelines/somatic_tumor_normal.php \
		-t_bam Sample_GS110168_03/GS110168_03.cram -n_bam Sample_GS110169_03/GS110169_03.cram \
		-out_folder out_case2 \
		-prefix GS110168_03-GS110169_03 \
		-steps vc,an,msi,an_rna\
		-t_rna_bam Sample_RX000001_01/RX000001_01.bam \
		-rna_ref_tissue gallbladder \
		-skip_correlation \
		-system system_pair.ini -n_system system_pair.ini \
		-min_depth_t 100 \
		-min_depth_n 100 \
		--log out_case2/GS110168_03-GS110169_03.log

tn_diff:
	# small variants
	php ../compare_variants.php -in1 expected_case2/GS110168_03-GS110169_03.GSvar -in2 out_case2/GS110168_03-GS110169_03.GSvar -add tumor_af,tumor_dp,normal_af,normal_dp,rna_ref_tpm
	# low coverge file
	php ../compare_bed.php expected_case2/GS110168_03-GS110169_03_stat_lowcov.bed out_case2/GS110168_03-GS110169_03_stat_lowcov.bed
	# qcML files
	php ../compare_qcml_xml.php expected_case2/GS110168_03-GS110169_03_stats_som.qcML out_case2/GS110168_03-GS110169_03_stats_som.qcML
	php ../compare_qcml_xml.php expected_case2/GS110168_03-GS110169_03_stats_other.qcML out_case2/GS110168_03-GS110169_03_stats_other.qcML
	# MSI file
	diff expected_case2/GS110168_03-GS110169_03_msi.tsv out_case2/GS110168_03-GS110169_03_msi.tsv
	# bafs
	diff expected_case2/GS110168_03-GS110169_03_bafs.igv out_case2/GS110168_03-GS110169_03_bafs.igv
	# cfDNA candidates
	diff --ignore-matching-lines '^##' expected_case2/GS110168_03-GS110169_03_cfDNA_candidates_monitoring.vcf out_case2/GS110168_03-GS110169_03_cfDNA_candidates/monitoring.vcf
	diff expected_case2/GS110168_03-GS110169_03_cfDNA_candidates_monitoring.tsv out_case2/GS110168_03-GS110169_03_cfDNA_candidates/monitoring.tsv

#Test MANTA SV calls with positive fusion control
tumor_normal_with_fusion: tn2_somatic tn2_diff

tn2_somatic:
	php ../../src/Pipelines/somatic_tumor_normal.php \
		-t_bam ../data/manta_sv_tumor.bam -n_bam ../data/manta_sv_normal.bam \
		-out_folder out_case3 \
		-prefix manta_sv_tumor-manta_sv_normal \
		-steps vc\
		-system system_pair_sv.ini -n_system system_pair_sv.ini \
		-skip_correlation \
		-min_depth_t 100 \
		-min_depth_n 100 \
		--log out_case3/tumor-normal.log
	rm ../data/manta_sv_normal_hla_genotyper.tsv ../data/manta_sv_tumor_hla_genotyper.tsv

tn2_diff:
	# remove gene annotation if present
	cut -f1-27 out_case3/manta_sv_tumor-manta_sv_normal_manta_var_structural.bedpe > out_case3/manta_sv_tumor-manta_sv_normal_manta_var_structural_cut.bedpe
	diff expected_case3/manta_sv_tumor-manta_sv_normal_manta_var_structural.bedpe out_case3/manta_sv_tumor-manta_sv_normal_manta_var_structural_cut.bedpe --ignore-matching-lines="##cmdline" --ignore-matching-lines="##"

create_test_case1:
	samtools view -bh /mnt/storage2/projects/diagnostic/SomaticAndTreatment/Sample_DX205723_01/DX205723_01.bam chr7:140723402-140929330 > test_case1_tmp.bam
	BamToFastq -in test_case1_tmp.bam -out1 som_pipeline_test_case1_L001_R1_001.fastq.gz -out2 som_pipeline_test_case1_L001_R2_001.fastq.gz
	rm test_case1_tmp.bam
	BedIntersect -in2 /mnt/storage2/megSAP/data/enrichment/ssSC_v5.bed  <<< "chr7	140723402	140929330" | BedSort | BedMerge > som_pipeline_test_case1_target.bed

#Test DeepSomatic as somatic variant caller	
tn_deepsom:
	mkdir -p out_case_deepsom
	rm -rf out_case_deepsom/GS110168_03-GS110169_03.log
	php ../../src/Pipelines/somatic_tumor_normal.php \
		-t_bam Sample_GS110168_03/GS110168_03.cram -n_bam Sample_GS110169_03/GS110169_03.cram \
		-out_folder out_case_deepsom \
		-prefix GS110168_03-GS110169_03 \
		-steps vc,an,msi,an_rna\
		-t_rna_bam Sample_RX000001_01/RX000001_01.bam \
		-rna_ref_tissue gallbladder \
		-skip_correlation \
		-system system_pair.ini -n_system system_pair.ini \
		-min_depth_t 100 \
		-min_depth_n 100 \
		-use_deepsomatic_test \
		--log out_case_deepsom/GS110168_03-GS110169_03.log
		
		
tn_deepsom_diff:
	# small variants
	php ../compare_variants.php -in1 expected_case_deepsom/GS110168_03-GS110169_03.GSvar -in2 out_case_deepsom/GS110168_03-GS110169_03.GSvar -add tumor_af,tumor_dp,normal_af,normal_dp,rna_ref_tpm
	# low coverge file
	php ../compare_bed.php expected_case_deepsom/GS110168_03-GS110169_03_stat_lowcov.bed out_case_deepsom/GS110168_03-GS110169_03_stat_lowcov.bed
	# qcML files
	php ../compare_qcml_xml.php expected_case_deepsom/GS110168_03-GS110169_03_stats_som.qcML out_case_deepsom/GS110168_03-GS110169_03_stats_som.qcML
	php ../compare_qcml_xml.php expected_case_deepsom/GS110168_03-GS110169_03_stats_other.qcML out_case_deepsom/GS110168_03-GS110169_03_stats_other.qcML
	# MSI file
	diff expected_case_deepsom/GS110168_03-GS110169_03_msi.tsv out_case_deepsom/GS110168_03-GS110169_03_msi.tsv
	# bafs
	diff expected_case_deepsom/GS110168_03-GS110169_03_bafs.igv out_case_deepsom/GS110168_03-GS110169_03_bafs.igv
	# cfDNA candidates
	diff --ignore-matching-lines '^##' expected_case_deepsom/GS110168_03-GS110169_03_cfDNA_candidates_monitoring.vcf out_case_deepsom/GS110168_03-GS110169_03_cfDNA_candidates/monitoring.vcf
	diff expected_case_deepsom/GS110168_03-GS110169_03_cfDNA_candidates_monitoring.tsv out_case_deepsom/GS110168_03-GS110169_03_cfDNA_candidates/monitoring.tsv