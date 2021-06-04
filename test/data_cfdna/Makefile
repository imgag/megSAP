all: clean personalized cfdna

clean:
	rm -rf CFDNA_01/ CFDNA_02/

personalized: link_p map_p vc_p an_p diff_p

link_p:
	mkdir -p CFDNA_01
	php ../link_test_data.php FO13609x01_01_L001_R1_001.fastq.gz CFDNA_01/CFDNA_01_L001_R1_001.fastq.gz
	php ../link_test_data.php FO13609x01_01_L001_R2_001.fastq.gz CFDNA_01/CFDNA_01_L001_R2_001.fastq.gz

map_p:
	php ../../src/Pipelines/analyze_cfdna.php -folder CFDNA_01 -name CFDNA_01 -system CFDNA_01.ini -steps ma -threads 1 -tumor_id TumorID -roi_patient roi_patient.bed

vc_p:
	php ../../src/Pipelines/analyze_cfdna.php -folder CFDNA_01 -name CFDNA_01 -system CFDNA_01.ini -steps vc -tumor_id TumorID -roi_patient roi_patient.bed

an_p:
	php ../../src/Pipelines/analyze_cfdna.php -folder CFDNA_01 -name CFDNA_01 -system CFDNA_01.ini -steps an -threads 1 -tumor_id TumorID -roi_patient roi_patient.bed

diff_p:
	php ../compare_qcml_xml.php CFDNA_01_stats_map.qcML.out CFDNA_01/CFDNA_01_stats_map.qcML 0.02
	diff --ignore-matching-lines '^##' CFDNA_01_var.vcf.out CFDNA_01/CFDNA_01_var.vcf
	php ../compare_variants.php -in1 CFDNA_01.GSvar.out -in2 CFDNA_01/CFDNA_01.GSvar


cfdna:  link_c map_c vc_c an_c diff_c

link_c:
	mkdir -p CFDNA_02
	php ../link_test_data.php FO13609x01_01_L001_R1_001.fastq.gz CFDNA_02/CFDNA_02_L001_R1_001.fastq.gz
	php ../link_test_data.php FO13609x01_01_L001_R2_001.fastq.gz CFDNA_02/CFDNA_02_L001_R2_001.fastq.gz

map_c:
	php ../../src/Pipelines/analyze_cfdna.php -folder CFDNA_02 -name CFDNA_02 -system CFDNA_02.ini -steps ma -threads 1

vc_c:
	php ../../src/Pipelines/analyze_cfdna.php -folder CFDNA_02 -name CFDNA_02 -system CFDNA_02.ini -steps vc

an_c:
	php ../../src/Pipelines/analyze_cfdna.php -folder CFDNA_02 -name CFDNA_02 -system CFDNA_02.ini -steps an

diff_c:
	php ../compare_qcml_xml.php CFDNA_02_stats_map.qcML.out CFDNA_02/CFDNA_02_stats_map.qcML 0.02
	diff --ignore-matching-lines '^##' CFDNA_02_var.vcf.out CFDNA_02/CFDNA_02_var.vcf
	php ../compare_variants.php -in1 CFDNA_02.GSvar.out -in2 CFDNA_02/CFDNA_02.GSvar