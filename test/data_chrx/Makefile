
all: clean link map vc an cn diff

clean:
	rm -rf *.fastq.gz *.log *.txt *.GSvar *.vcf.gz *.vcf.gz.tbi *.bam *.bai *.qcML *.bed *_cnvs.seg *_cnvs.tsv

link:
	php ../link_test_data.php NA12878_17_GTAGAGGATATCCTCT_L001_R1_001.fastq.gz
	php ../link_test_data.php NA12878_17_GTAGAGGATATCCTCT_L001_R2_001.fastq.gz

map:
	rm -rf NA12878_17_log1_map.log
	php ../../src/Pipelines/analyze.php  -folder . -name NA12878_17  -system system.ini -steps ma -threads 1

vc:
	rm -rf NA12878_17_log2_vc.log
	php ../../src/Pipelines/analyze.php  -folder . -name NA12878_17  -system system.ini -steps vc

an:
	rm -rf NA12878_17_log3_anno.log
	rm -rf NA12878_17_log5_diag.log
	php ../../src/Pipelines/analyze.php  -folder . -name NA12878_17  -system system.ini -steps an

cn:
	php ../../src/Pipelines/analyze.php -folder . -name NA12878_17 -system system.ini -steps cn

diff:
	php ../compare_variants.php -in1 NA12878_17.GSvar.out -in2 NA12878_17.GSvar
	php ../compare_qcml_xml.php NA12878_17_stats_fastq.qcML.out NA12878_17_stats_fastq.qcML
	php ../compare_qcml_xml.php NA12878_17_stats_map.qcML.out NA12878_17_stats_map.qcML 0.01
	php ../compare_qcml_xml.php NA12878_17_stats_vc.qcML.out NA12878_17_stats_vc.qcML
	diff NA12878_17_cnvs.tsv.out NA12878_17_cnvs.tsv
	numdiff NA12878_17_cnvs.seg.out NA12878_17_cnvs.seg --absolute-tolerance 0.04
	php ../compare_bed.php NA12878_17_nxLRRK2_SNCA_lowcov.bed.out NA12878_17_nxLRRK2_SNCA_lowcov.bed


