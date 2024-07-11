all: clean personalized cfdna

clean:
	rm -rf CFDNA_01/ CFDNA_02/

personalized: link_p map_p vc_p diff_p

link_p:
	mkdir -p CFDNA_01
	php ../link_test_data.php FO13609x01_01_L001_R1_001.fastq.gz CFDNA_01/CFDNA_01_L001_R1_001.fastq.gz
	php ../link_test_data.php FO13609x01_01_L001_R2_001.fastq.gz CFDNA_01/CFDNA_01_L001_R2_001.fastq.gz

map_p:
	php ../../src/Pipelines/analyze_cfdna.php -folder CFDNA_01 -name CFDNA_01 -system CFDNA_01.ini -steps ma -skip_tumor -target roi_patient.bed -monitoring_vcf CFDNA_01_monitoring.vcf -threads 1

vc_p:
	php ../../src/Pipelines/analyze_cfdna.php -folder CFDNA_01 -name CFDNA_01 -system CFDNA_01.ini -steps vc -skip_tumor -target roi_patient.bed -monitoring_vcf CFDNA_01_monitoring.vcf -threads 1

diff_p:
	php ../compare_qcml_xml.php CFDNA_01_stats_map.qcML.out CFDNA_01/CFDNA_01_stats_map.qcML 0.02
	php ../compare_qcml_xml.php CFDNA_01_stats_fastq.qcML.out CFDNA_01/CFDNA_01_stats_fastq.qcML 0.02
	php ../compare_qcml_xml.php CFDNA_01_stats_cfDNA.qcML.out CFDNA_01/CFDNA_01_stats_cfDNA.qcML 0.02
	diff --ignore-matching-lines '^##' CFDNA_01_var.vcf.out CFDNA_01/CFDNA_01_var.vcf
	php ../compare_variants.php -in1 CFDNA_01.GSvar.out -in2 CFDNA_01/CFDNA_01.GSvar
	diff CFDNA_01_bg.mrd.out CFDNA_01/umiVar/CFDNA_01_bg.mrd
	diff CFDNA_01_bg_unfiltered.mrd.out CFDNA_01/umiVar/CFDNA_01_bg_unfiltered.mrd
	diff CFDNA_01_monitoring_counts.tsv.out CFDNA_01/umiVar/CFDNA_01_monitoring_counts.tsv
	diff CFDNA_01_monitoring_counts_unfiltered.tsv.out CFDNA_01/umiVar/CFDNA_01_monitoring_counts_unfiltered.tsv

ra_p:
	mv CFDNA_01/CFDNA_01.GSvar CFDNA_01/CFDNA_01_beforeRA.GSvar
	php ../../src/Pipelines/analyze_cfdna.php -folder CFDNA_01 -name CFDNA_01 -system CFDNA_01.ini -steps vc -skip_tumor -target roi_patient.bed -monitoring_vcf CFDNA_01_monitoring.vcf -threads 1 -annotation_only
	diff CFDNA_01/CFDNA_01.GSvar CFDNA_01/CFDNA_01_beforeRA.GSvar

cfdna:  link_c map_c vc_c diff_c

link_c:
	mkdir -p CFDNA_02
	php ../link_test_data.php FO13609x01_01_L001_R1_001.fastq.gz CFDNA_02/CFDNA_02_L001_R1_001.fastq.gz
	php ../link_test_data.php FO13609x01_01_L001_R2_001.fastq.gz CFDNA_02/CFDNA_02_L001_R2_001.fastq.gz

map_c:
	php ../../src/Pipelines/analyze_cfdna.php -folder CFDNA_02 -name CFDNA_02 -system CFDNA_02.ini -skip_tumor -steps ma -threads 1

vc_c:
	php ../../src/Pipelines/analyze_cfdna.php -folder CFDNA_02 -name CFDNA_02 -system CFDNA_02.ini -skip_tumor -steps vc -threads 1

diff_c:
	php ../compare_qcml_xml.php CFDNA_02_stats_map.qcML.out CFDNA_02/CFDNA_02_stats_map.qcML 0.02
	php ../compare_qcml_xml.php CFDNA_02_stats_fastq.qcML.out CFDNA_02/CFDNA_02_stats_fastq.qcML 0.02
	php ../compare_qcml_xml.php CFDNA_02_stats_cfDNA.qcML.out CFDNA_02/CFDNA_02_stats_cfDNA.qcML 0.02
	diff --ignore-matching-lines '^##' CFDNA_02_var.vcf.out CFDNA_02/CFDNA_02_var.vcf
	php ../compare_variants.php -in1 CFDNA_02.GSvar.out -in2 CFDNA_02/CFDNA_02.GSvar

ra_c:
	mv CFDNA_02/CFDNA_02.GSvar CFDNA_02/CFDNA_02_beforeRA.GSvar
	php ../../src/Pipelines/analyze_cfdna.php -folder CFDNA_02 -name CFDNA_02 -system CFDNA_02.ini -steps vc -threads 1 -annotation_only
	diff CFDNA_02/CFDNA_02.GSvar CFDNA_02/CFDNA_02_beforeRA.GSvar