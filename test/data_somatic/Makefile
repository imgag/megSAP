#TODO
all: clear link somatic_single somatic_pair somatic_amp somatic_uid diff_single diff_pair diff_amp diff_uid

clear:
	rm -rf *.bam *.bai *.log *.txt *.GSvar *.vcf *qci.vcf.gz *.fastq.gz *.qcML
	rm -rf Sample_GS110168
	rm -rf Sample_GS110169
	rm -rf Sample_GS150344_01
	rm -rf Sample_SampleUID1_01
	rm -rf Sample_SampleUID2_01
	rm -rf results_single
	rm -rf results_pair1
	rm -rf results_pair2
	rm -rf results_uid1
	rm -rf results_uid2
	rm -rf results_amp

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
	mkdir Sample_SampleUID1_01
	php ../link_test_data.php SampleUID1_01_TTGACA_L004_R1_001.fastq.gz Sample_SampleUID1_01/SampleUID1_01_TTGACA_L004_R1_001.fastq.gz
	php ../link_test_data.php SampleUID1_01_TTGACA_L004_R2_001.fastq.gz Sample_SampleUID1_01/SampleUID1_01_TTGACA_L004_R2_001.fastq.gz
	mkdir Sample_SampleUID2_01
	php ../link_test_data.php SampleUID2_01_GATAGA_L001_index_001.fastq.gz Sample_SampleUID2_01/SampleUID2_01_GATAGA_L001_index_001.fastq.gz
	php ../link_test_data.php SampleUID2_01_GATAGA_L001_R1_001.fastq.gz Sample_SampleUID2_01/SampleUID2_01_GATAGA_L001_R1_001.fastq.gz
	php ../link_test_data.php SampleUID2_01_GATAGA_L001_R2_001.fastq.gz Sample_SampleUID2_01/SampleUID2_01_GATAGA_L001_R2_001.fastq.gz
	mkdir results_amp
	mkdir results_uid1
	mkdir results_uid2

somatic_single:
	php ../../src/Pipelines/somatic_dna.php -p_folder . -t_id GS150344_01 -n_id na -o_folder results_single -steps ma,vc,an -t_sys system_single.ini --log results_single/GS150344.log

diff_single:
	php ../compare_variants.php GS150344_01.GSvar.out_freebayes results_single/GS150344_01.GSvar

somatic_pair:
	#strelka
	php ../../src/Pipelines/somatic_dna.php -p_folder . -t_id GS110168 -n_id GS110169 -o_folder results_pair1 -steps ma,vc,an,ci -reduce_variants_filter -t_sys system_pair.ini -n_sys system_pair.ini --log results_pair1/GS110168-GS110169.log
	#freebayes (previous mapping results are used!)
	php ../../src/Pipelines/somatic_dna.php -p_folder . -t_id GS110168 -n_id GS110169 -o_folder results_pair2 -steps vc,an -t_sys system_pair.ini -n_sys system_pair.ini -freebayes --log results_pair2/GS110168-GS110169.log

diff_pair:
	php ../compare_variants.php GS110168-GS110169.GSvar.out_strelka results_pair1/GS110168-GS110169.GSvar
	php ../compare_variants.php GS110168-GS110169.GSvar.out_freebayes results_pair2/GS110168-GS110169.GSvar
	ls Sample_*/*.qcML | sort | xargs grep "<qualityParameter" > sample_qcML.txt
	ls results_pair1/*.qcML | sort | xargs grep "<qualityParameter" > somatic_pair1_qcML.txt
	ls results_pair2/*.qcML | sort | xargs grep "<qualityParameter" > somatic_pair2_qcML.txt
	php ../compare_qcml.php sample_qcML.txt.out sample_qcML.txt
	php ../compare_qcml.php somatic_pair1_qcML.txt.out somatic_pair1_qcML.txt
	php ../compare_qcml.php somatic_pair2_qcML.txt.out somatic_pair2_qcML.txt
	php ../compare_variants.php GS110168-GS110169_var_qci.pair1.out.vcf.gz results_pair1/GS110168-GS110169_var_qci.vcf.gz

somatic_amp:
	#php ../../src/Pipelines/somatic_dna.php -p_folder . -t_id SampleUID2_01 -n_id na -o_folder results_amp -steps ma,vc,an -reduce_variants_filter -t_sys system_amp.ini --log results_amp/SampleUID2_01.log

diff_amp:
	#php ../compare_variants.php SampleUID2_01.GSvar.out_amp results_amp/SampleUID2_01.GSvar

somatic_uid:
	#MIPs
	#php ../../src/Pipelines/somatic_dna.php -p_folder . -t_id SampleUID1_01 -n_id na -o_folder results_uid1 -steps ma,vc,an -reduce_variants_filter -no_softclip -t_sys system_uid1.ini --log results_uid1/SampleUID1_01.log
	#HPHS
	#php ../../src/Pipelines/somatic_dna.php -p_folder . -t_id SampleUID2_01 -n_id na -o_folder results_uid2 -steps ma,vc,an -reduce_variants_filter -no_softclip -t_sys system_uid2.ini --log results_uid2/SampleUID2_01.log

diff_uid:
	#TODO: Comparison of SampleUID1.
	#php ../compare_variants.php SampleUID1_01.GSvar.out_uid1 results_uid1/SampleUID1_01.GSvar
	#TODO: Comparison of SampleUID2.
	#php ../compare_variants.php SampleUID2_01.GSvar.out_uid2 results_uid2/SampleUID2_01.GSvar
	
#debugging targets
anno1:
	php ../../src/Pipelines/analyze.php -folder . -name GS110168 -system system_pair.ini -steps an

anno2:
	php ../../src/Pipelines/analyze.php -folder . -name GS110169 -system system_pair.ini -steps an
