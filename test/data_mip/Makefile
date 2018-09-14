
all: clean link map vc an cn diff

clean:
	rm -rf *.fastq.gz *.log *.txt *.GSvar *.vcf.gz *.vcf.gz.tbi *.bam *.bai *.qcML *_cnvs.seg *_cnvs.tsv
	
link:
	php ../link_test_data.php NA12878_39_S51_L001_R1_001.fastq.gz
	php ../link_test_data.php NA12878_39_S51_L001_R2_001.fastq.gz
	php ../link_test_data.php NA12878_39_S51_L002_R1_001.fastq.gz
	php ../link_test_data.php NA12878_39_S51_L002_R2_001.fastq.gz
	php ../link_test_data.php NA12878_39_S51_L003_R1_001.fastq.gz
	php ../link_test_data.php NA12878_39_S51_L003_R2_001.fastq.gz
	php ../link_test_data.php NA12878_39_S51_L004_R1_001.fastq.gz
	php ../link_test_data.php NA12878_39_S51_L004_R2_001.fastq.gz


map:
	rm -rf NA12878_39_log1_map.log
	php ../../src/Pipelines/analyze.php -folder . -name NA12878_39 -system system.ini -steps ma -threads 1

vc:
	rm -rf NA12878_39_log2_vc.log
	php ../../src/Pipelines/analyze.php -folder . -name NA12878_39 -system system.ini -steps vc

an:
	rm -rf NA12878_39_log3_anno.log
	php ../../src/Pipelines/analyze.php -folder . -name NA12878_39 -system system.ini -steps an

cn:
	php ../../src/Pipelines/analyze.php -folder . -name NA12878_39 -system system.ini -steps cn

diff:
	php ../compare_variants.php -in1 NA12878_39.GSvar.out -in2 NA12878_39.GSvar