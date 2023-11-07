
all: clean link map vc cn diff reanno

clean:
	rm -rf *.fastq.gz *.log *.txt *.GSvar *.vcf.gz *.vcf.gz.tbi *.cram *.crai *.qcML *.bed *_cnvs_clincnv.seg *_cnvs_clincnv.tsv *.bedpe manta_evid *_ancestry.tsv *_bafs.igv
	
link:
	php ../link_test_data.php NA12878_02_AACGTGAT_L001_R1_001.fastq.gz
	php ../link_test_data.php NA12878_02_AACGTGAT_L001_R2_001.fastq.gz
	php ../link_test_data.php NA12878_03_AACGTGAT_L001_R1_001.fastq.gz
	php ../link_test_data.php NA12878_03_AACGTGAT_L001_R2_001.fastq.gz

map:
	php ../../src/Pipelines/analyze.php -folder . -name NA12878_03 -system system.ini -steps ma -threads 1

vc:
	php ../../src/Pipelines/analyze.php -folder . -name NA12878_03 -system system.ini -steps vc

cn:
	php ../../src/Pipelines/analyze.php -folder . -name NA12878_03 -system system.ini -steps cn

diff:
	php ../compare_variants.php -in1 NA12878_03.GSvar.out -in2 NA12878_03.GSvar
	php ../compare_qcml_xml.php NA12878_03_stats_fastq.qcML.out NA12878_03_stats_fastq.qcML
	php ../compare_qcml_xml.php NA12878_03_stats_map.qcML.out NA12878_03_stats_map.qcML
	php ../compare_qcml_xml.php NA12878_03_stats_vc.qcML.out NA12878_03_stats_vc.qcML
	php ../compare_qcml_xml.php NA12878_03_stats_other.qcML.out NA12878_03_stats_other.qcML
	diff NA12878_03_cnvs_clincnv.tsv.out NA12878_03_cnvs_clincnv.tsv --ignore-matching-lines="#"
	php ../compare_bed.php NA12878_03_hpHBOCv5_lowcov.bed.out NA12878_03_hpHBOCv5_lowcov.bed

reanno:
	cp NA12878_03.GSvar NA12878_03_old.GSvar
	cp NA12878_03_cnvs_clincnv.tsv NA12878_03_old_cnvs_clincnv.tsv
	php ../../src/Pipelines/analyze.php -folder . -name NA12878_03 -system system.ini -steps vc,cn -annotation_only
	diff NA12878_03.GSvar NA12878_03_old.GSvar
	diff NA12878_03_cnvs_clincnv.tsv NA12878_03_old_cnvs_clincnv.tsv
