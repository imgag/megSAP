
all: clean link map vc cn sv an diff

clean:
	rm -rf *.fastq.gz *.log *.txt *.GSvar *.vcf.gz *.vcf.gz.tbi *.bam *.bai *.qcML *.bed *_cnvs_clincnv.seg *_cnvs_clincnv.tsv *.bedpe clair_temp repeat_expansions *_bafs.igv *_circos.png *_mosaic_cnvs.tsv *_rohs.tsv *_ancestry.tsv

link:
	php ../link_test_data.php 23014LRa023L2_01_chr1_longread.fastq.gz

map:
	php ../../src/Pipelines/analyze_longread.php -folder . -name 23014LRa023L2_01 -system system.ini -steps ma -threads 1
	
vc:
	php ../../src/Pipelines/analyze_longread.php -folder . -name 23014LRa023L2_01 -system system.ini -steps vc 
	
cn:
	php ../../src/Pipelines/analyze_longread.php -folder . -name 23014LRa023L2_01 -system system.ini -steps cn

sv:
	php ../../src/Pipelines/analyze_longread.php -folder . -name 23014LRa023L2_01 -system system.ini -steps sv

an:
	php ../../src/Pipelines/analyze_longread.php -folder . -name 23014LRa023L2_01 -system system.ini -steps an
	# remove gene annotation
	cut -f1-22 23014LRa023L2_01_var_structural_variants.bedpe > 23014LRa023L2_01_var_structural_variants_cut.bedpe
	# create file from OMIM annotation and CNV overlap
	# cat 23014LRa023L2_01_var_structural_variants.bedpe | egrep -v "^##"| rev | cut -f 1,2| rev > 23014LRa023L2_01_var_structural_variants_cut2.bedpe

diff:
	php ../compare_variants.php -in1 23014LRa023L2_01.GSvar.out -in2 23014LRa023L2_01.GSvar
	php ../compare_qcml_xml.php 23014LRa023L2_01_stats_fastq.qcML.out 23014LRa023L2_01_stats_fastq.qcML
	php ../compare_qcml_xml.php 23014LRa023L2_01_stats_map.qcML.out 23014LRa023L2_01_stats_map.qcML 0.01
	php ../compare_qcml_xml.php 23014LRa023L2_01_stats_vc.qcML.out 23014LRa023L2_01_stats_vc.qcML
	php ../compare_qcml_xml.php 23014LRa023L2_01_stats_other.qcML.out 23014LRa023L2_01_stats_other.qcML
	diff 23014LRa023L2_01_cnvs_clincnv.tsv.out 23014LRa023L2_01_cnvs_clincnv.tsv --ignore-matching-lines="#"
	# numdiff 23014LRa023L2_01_cnvs_clincnv.seg.out 23014LRa023L2_01_cnvs_clincnv.seg --absolute-tolerance 0.04
	# php ../compare_bed.php 23014LRa023L2_01_nxLRRK2_SNCA_lowcov.bed.out 23014LRa023L2_01_nxLRRK2_SNCA_lowcov.bed
	diff 23014LRa023L2_01_var_structural_variants.bedpe.out 23014LRa023L2_01_var_structural_variants_cut.bedpe --ignore-matching-lines="^##"
	# diff 23014LRa023L2_01_var_structural_variants_last_cols.bedpe.out 23014LRa023L2_01_var_structural_variants_cut2.bedpe --ignore-matching-lines="^##"
	diff 23014LRa023L2_01_prs.tsv.out 23014LRa023L2_01_prs.tsv

