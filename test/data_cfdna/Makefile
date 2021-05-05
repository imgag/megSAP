all: clean link map vc an diff

clean:
	rm -rf CFDNA_01*.fastq.gz CFDNA_01*.bam* CFDNA_01*.qcML CFDNA_01*.bed CFDNA_01*_var.vcf CFDNA_01*_var.tsv CFDNA_01*_var.mrd *.log *.GSvar umiVar_tmp*

link:
	php ../link_test_data.php FO13609x01_01_L001_R1_001.fastq.gz CFDNA_01_L001_R1_001.fastq.gz
	php ../link_test_data.php FO13609x01_01_L001_R2_001.fastq.gz CFDNA_01_L001_R2_001.fastq.gz

map:
	php ../../src/Pipelines/analyze_cfdna.php -folder . -name CFDNA_01 -system system.ini -steps ma -threads 1 -tumor_id TumorID -roi_patient roi_patient.bed

vc:
	php ../../src/Pipelines/analyze_cfdna.php -folder . -name CFDNA_01 -system system.ini -steps vc -tumor_id TumorID -roi_patient roi_patient.bed

an:
	php ../../src/Pipelines/analyze_cfdna.php -folder . -name CFDNA_01 -system system.ini -steps an -threads 1 -tumor_id TumorID -roi_patient roi_patient.bed

diff:
	php ../compare_qcml_xml.php CFDNA_01_stats_map.qcML.out CFDNA_01_stats_map.qcML 0.02
	diff --ignore-matching-lines '^##' CFDNA_01_var.vcf.out CFDNA_01_var.vcf
	php ../compare_variants.php -in1 CFDNA_01.GSvar.out -in2 CFDNA_01.GSvar
