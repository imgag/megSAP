
all: clean link map vc an cn diff

clean:
	rm -rf *.fastq.gz *.log *.txt *.GSvar *.vcf.gz *.vcf.gz.tbi *.bam *.bai *.qcML *.bed *_cnvs.seg *_cnvs.tsv
	
link:
	php ../link_test_data.php NA12878_02_AACGTGAT_L001_R1_001.fastq.gz
	php ../link_test_data.php NA12878_02_AACGTGAT_L001_R2_001.fastq.gz
	php ../link_test_data.php NA12878_03_AACGTGAT_L001_R1_001.fastq.gz
	php ../link_test_data.php NA12878_03_AACGTGAT_L001_R2_001.fastq.gz


map:
	rm -rf NA12878_03_log1_map.log
	php ../../src/Pipelines/analyze.php -folder . -name NA12878_03 -system system.ini -steps ma -threads 1

vc:
	rm -rf NA12878_03_log2_vc.log
	php ../../src/Pipelines/analyze.php -folder . -name NA12878_03 -system system.ini -steps vc

an:
	rm -rf NA12878_03_log3_anno.log
	php ../../src/Pipelines/analyze.php -folder . -name NA12878_03 -system system.ini -steps an

cn:
	php ../../src/Pipelines/analyze.php -folder . -name NA12878_03 -system system.ini -steps cn

diff:
	php ../compare_variants.php -in1 NA12878_03.GSvar.out -in2 NA12878_03.GSvar
	php ../compare_qcml_xml.php NA12878_03_stats_fastq.qcML.out NA12878_03_stats_fastq.qcML
	php ../compare_qcml_xml.php NA12878_03_stats_map.qcML.out NA12878_03_stats_map.qcML
	php ../compare_qcml_xml.php NA12878_03_stats_vc.qcML.out NA12878_03_stats_vc.qcML
	diff NA12878_03_cnvs.tsv.out NA12878_03_cnvs.tsv --ignore-matching-lines="CnvHunter version:"
	php ../compare_bed.php NA12878_03_hpHBOCv5_lowcov.bed.out NA12878_03_hpHBOCv5_lowcov.bed
