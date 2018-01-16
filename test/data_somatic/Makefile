all: clear link somatic_single1 somatic_pair1 somatic_pair2

clear:
	rm -rf *.bam *.bai *.log *.txt *.GSvar *.vcf *qci.vcf.gz *.fastq.gz *.qcML
	rm -rf Sample_GS110168 Sample_GS110169 Sample_GS150344_01 Sample_Tumor2 Sample_Normal2 results_single results_pair1 results_pair2 results_amp Somatic_Tumor2-Normal2

link:
	mkdir -p Sample_GS110169 Sample_GS110168 Sample_GS150344_01 Sample_GS150344_01 Sample_Tumor2 Sample_Normal2 Somatic_Tumor2-Normal2 results_single results_pair1 results_pair2
	php ../link_test_data.php GS110168_chr21_L007_R1_001.fastq.gz Sample_GS110168/GS110168_chr21_L007_R1_001.fastq.gz
	php ../link_test_data.php GS110168_chr21_L007_R2_001.fastq.gz Sample_GS110168/GS110168_chr21_L007_R2_001.fastq.gz
	php ../link_test_data.php GS110169_chr21_L008_R1_001.fastq.gz Sample_GS110169/GS110169_chr21_L008_R1_001.fastq.gz
	php ../link_test_data.php GS110169_chr21_L008_R2_001.fastq.gz Sample_GS110169/GS110169_chr21_L008_R2_001.fastq.gz
	php ../link_test_data.php GS150344_01_GCTCGGTAextracted_L001_R1_001.fastq.gz Sample_GS150344_01/GS150344_01_GCTCGGTAextracted_L001_R1_001.fastq.gz
	php ../link_test_data.php GS150344_01_GCTCGGTAextracted_L001_R2_001.fastq.gz Sample_GS150344_01/GS150344_01_GCTCGGTAextracted_L001_R2_001.fastq.gz
	php ../link_test_data.php GS150344_01_GCTCGGTAextracted_L001_R1_001.fastq.gz Sample_GS150344_01/GS150344_01_GCTCGGTAextracted_L001_R1_001.fastq.gz
	php ../link_test_data.php GS150344_01_GCTCGGTAextracted_L001_R2_001.fastq.gz Sample_GS150344_01/GS150344_01_GCTCGGTAextracted_L001_R2_001.fastq.gz
	php ../link_test_data.php Tumor2_XX_L00X_R1_001.fastq.gz Sample_Tumor2/Tumor2_XX_L00X_R1_001.fastq.gz
	php ../link_test_data.php Tumor2_XX_L00X_R2_001.fastq.gz Sample_Tumor2/Tumor2_XX_L00X_R2_001.fastq.gz
	php ../link_test_data.php Normal2_XX_L00X_R1_001.fastq.gz Sample_Normal2/Normal2_XX_L00X_R1_001.fastq.gz
	php ../link_test_data.php Normal2_XX_L00X_R2_001.fastq.gz Sample_Normal2/Normal2_XX_L00X_R2_001.fastq.gz

somatic_single1:
	#single1
	php ../../src/Pipelines/somatic_dna.php -p_folder . -t_id GS150344_01 -n_id na -o_folder results_single -steps ma,vc,an -t_sys system_single.ini --log results_single/GS150344.log
	php ../compare_variants.php GS150344_01.GSvar.out_freebayes results_single/GS150344_01.GSvar

somatic_pair1:
	#pair1
	#strelka
	php ../../src/Pipelines/somatic_dna.php -p_folder . -t_id GS110168 -n_id GS110169 -o_folder results_pair1 -steps ma,vc,an,ci -reduce_variants_filter -t_sys system_pair.ini -n_sys system_pair.ini --log results_pair1/GS110168-GS110169.log
	php ../compare_variants.php GS110168-GS110169.GSvar.out_strelka results_pair1/GS110168-GS110169.GSvar
	ls Sample_*/*.qcML | sort | xargs grep "<qualityParameter" > sample_qcML.txt
	ls results_pair1/*.qcML | sort | xargs grep "<qualityParameter" > somatic_pair1_qcML.txt
	php ../compare_qcml.php sample_qcML.txt.out sample_qcML.txt
	php ../compare_qcml.php somatic_pair1_qcML.txt.out somatic_pair1_qcML.txt
	php ../compare_variants.php GS110168-GS110169_var_qci.pair1.out.vcf.gz results_pair1/GS110168-GS110169_var_qci.vcf.gz
	test -s results_pair1/GS110168-GS110169_cgi_drug_prescription.tsv
	test -s results_pair1/GS110168-GS110169_cgi_drug_prescription_bioactivities.tsv
	test -s results_pair1/GS110168-GS110169_cgi_mutation_analysis.tsv
	#pair1
	#freebayes (previous mapping results are used!)
	php ../../src/Pipelines/somatic_dna.php -p_folder . -t_id GS110168 -n_id GS110169 -o_folder results_pair2 -steps vc,an,ci -t_sys system_pair.ini -n_sys system_pair.ini -freebayes --log results_pair2/GS110168-GS110169.log
	ls Sample_*/*.qcML | sort | xargs grep "<qualityParameter" > sample_qcML.txt
	ls results_pair2/*.qcML | sort | xargs grep "<qualityParameter" > somatic_pair2_qcML.txt
	php ../compare_variants.php GS110168-GS110169.GSvar.out_freebayes results_pair2/GS110168-GS110169.GSvar
	php ../compare_qcml.php sample_qcML.txt.out sample_qcML.txt
	php ../compare_qcml.php somatic_pair2_qcML.txt.out somatic_pair2_qcML.txt
	test -s results_pair2/GS110168-GS110169_cgi_drug_prescription.tsv
	test -s results_pair2/GS110168-GS110169_cgi_drug_prescription_bioactivities.tsv
	test -s results_pair2/GS110168-GS110169_cgi_mutation_analysis.tsv

somatic_pair2:
	#pair2
	#strelka
	php ../../src/Pipelines/somatic_dna.php -p_folder . -t_id Tumor2 -n_id Normal2 -o_folder Somatic_Tumor2-Normal2 -steps ma,vc,an,ci -t_sys system_pair.ini -n_sys system_pair.ini --log Somatic_Tumor2-Normal2/Tumor2-Normal2.log
	php ../compare_variants.php Tumor2-Normal2_var_qci.out.vcf.gz Somatic_Tumor2-Normal2/Tumor2-Normal2_var_qci.vcf.gz
	php ../compare_variants.php Tumor2-Normal2.GSvar.out Somatic_Tumor2-Normal2/Tumor2-Normal2.GSvar
	ls Sample_Tumor2/*.qcML | sort | xargs grep "<qualityParameter" > somatic_Tumor2-Normal2_qcML.txt
	ls Sample_Normal2/*.qcML | sort | xargs grep "<qualityParameter" > somatic_Tumor2-Normal2_qcML.txt
	ls Somatic_Tumor2-Normal2/*.qcML | sort | xargs grep "<qualityParameter" > somatic_Tumor2-Normal2_qcML.txt
	php ../compare_qcml.php somatic_Tumor2-Normal2_qcML.txt.out somatic_Tumor2-Normal2_qcML.txt
	test -s Somatic_Tumor2-Normal2/Tumor2-Normal2_cgi_drug_prescription.tsv
	test -s Somatic_Tumor2-Normal2/Tumor2-Normal2_cgi_drug_prescription_bioactivities.tsv
	test -s Somatic_Tumor2-Normal2/Tumor2-Normal2_cgi_mutation_analysis.tsv