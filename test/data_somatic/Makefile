
all: clear link somatic_single somatic_pair diff_single diff_pair 

clear:
	rm -rf *.bam *.bai *.log *.txt *.GSvar *.vcf *.fastq.gz *.qcML
	rm -rf Sample_GS110168
	rm -rf Sample_GS110169
	rm -rf Sample_GS150344_01
	rm -rf results_single
	rm -rf results_pair1
	rm -rf results_pair2

link:
	mkdir -p Sample_GS110168
	php ../link_test_data.php GS110168_chr21_L007_R1_001.fastq.gz Sample_GS110168/GS110168_chr21_L007_R1_001.fastq.gz
	php ../link_test_data.php GS110168_chr21_L007_R2_001.fastq.gz Sample_GS110168/GS110168_chr21_L007_R2_001.fastq.gz
	mkdir -p Sample_GS110169
	php ../link_test_data.php GS110169_chr21_L008_R1_001.fastq.gz Sample_GS110169/GS110169_chr21_L008_R1_001.fastq.gz
	php ../link_test_data.php GS110169_chr21_L008_R2_001.fastq.gz Sample_GS110169/GS110169_chr21_L008_R2_001.fastq.gz
	mkdir -p Sample_GS150344_01
	php ../link_test_data.php GS150344_01_GCTCGGTAextracted_L001_R1_001.fastq.gz Sample_GS150344_01/GS150344_01_GCTCGGTAextracted_L001_R1_001.fastq.gz
	php ../link_test_data.php GS150344_01_GCTCGGTAextracted_L001_R2_001.fastq.gz Sample_GS150344_01/GS150344_01_GCTCGGTAextracted_L001_R2_001.fastq.gz
	mkdir results_single
	mkdir results_pair1
	mkdir results_pair2

somatic_single:
	php ../../src/Pipelines/somatic_dna.php -p_folder . -t_id GS150344_01 -n_id na -o_folder results_single -steps ma,vc,an -t_sys system_single.ini --log results_single/GS150344.log

diff_single:
	php ../compare_variants.php GS150344_01.GSvar.out_freebayes results_single/GS150344_01.GSvar

somatic_pair:
	#strelka
	php ../../src/Pipelines/somatic_dna.php -p_folder . -t_id GS110168 -n_id GS110169 -o_folder results_pair1 -steps ma,vc,an -reduce_variants_filter -t_sys system_pair.ini -n_sys system_pair.ini --log results_pair1/GS110168-GS110169.log
	#freebayes (previous mapping results are used!)
	php ../../src/Pipelines/somatic_dna.php -p_folder . -t_id GS110168 -n_id GS110169 -o_folder results_pair2 -steps vc,an -t_sys system_pair.ini -n_sys system_pair.ini -freebayes --log results_pair2/GS110168-GS110169.log

diff_pair:
	php ../compare_variants.php GS110168-GS110169.GSvar.out_strelka results_pair1/GS110168-GS110169.GSvar
	php ../compare_variants.php GS110168-GS110169.GSvar.out_freebayes results_pair2/GS110168-GS110169.GSvar

#debugging targets
anno1:
	php ../../src/Pipelines/analyze.php -folder . -name GS110168 -system system_pair.ini -steps an

anno2:
	php ../../src/Pipelines/analyze.php -folder . -name GS110169 -system system_pair.ini -steps an
