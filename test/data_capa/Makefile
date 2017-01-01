
all: clear link
	make somatic
	make diff

clear:
	rm -rf *.bam *.bai *.log *.txt *.GSvar *.vcf *.fastq.gz *.qcML
	rm -rf Sample_GS150496_01
	rm -rf Sample_GS150497_01
	rm -rf Somatic_GS150497_01-GS150496_01

link:
	mkdir -p Sample_GS150496_01
	php ../link_test_data.php GS150496_01_S1_L001extracted_R1_001.fastq.gz Sample_GS150496_01/GS150496_01_S1_L001extracted_R1_001.fastq.gz
	php ../link_test_data.php GS150496_01_S1_L001extracted_R2_001.fastq.gz Sample_GS150496_01/GS150496_01_S1_L001extracted_R2_001.fastq.gz
	mkdir -p Sample_GS150497_01
	php ../link_test_data.php GS150497_01_S2_L001extracted_R1_001.fastq.gz Sample_GS150497_01/GS150497_01_S2_L001extracted_R1_001.fastq.gz
	php ../link_test_data.php GS150497_01_S2_L001extracted_R2_001.fastq.gz Sample_GS150497_01/GS150497_01_S2_L001extracted_R2_001.fastq.gz

somatic:
	mkdir Somatic_GS150497_01-GS150496_01
	php ../../src/Pipelines/somatic_capa.php -t_id GS150497_01 -t_sys system.ini -n_id GS150496_01 -n_sys system.ini -p_folder . -o_folder Somatic_GS150497_01-GS150496_01 -no_db -nsc --log GS150497_01-GS150496_01.log
	cp Somatic_GS150497_01-GS150496_01/GS150497_01-GS150496_01.GSvar GS150497_01-GS150496_01.GSvar
	cp Somatic_GS150497_01-GS150496_01/GS150497_01-GS150496_01_report.txt GS150497_01-GS150496_01_report.txt

diff:
	php ../compare_variants.php GS150497_01-GS150496_01.GSvar.out GS150497_01-GS150496_01.GSvar
	diff -b  -I '^Datum' -I '^Revision' -I '^Tumor:' -I '^Normal:' -I '^  Target' GS150497_01-GS150496_01_report.txt.out GS150497_01-GS150496_01_report.txt
