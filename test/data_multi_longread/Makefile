all: clean link call diff

clean:
	rm -rf *.fastq.gz *.log *.txt *.GSvar *.vcf.gz *.vcf.gz.tbi *.cram *.crai *.bam *.bai *.qcML *.bed *_cnvs_clincnv.seg *_cnvs_clincnv.tsv *.bedpe clair_temp repeat_expansions *_bafs.igv *_circos.png *_mosaic_cnvs.tsv *_rohs.tsv *_ancestry.tsv *.gvcf.gz *.gvcf.gz.tbi

link:
	php ../link_test_data.php 24067LRa008_01.cram
	php ../link_test_data.php 24067LRa008_01.cram.crai
	php ../link_test_data.php 24067LRa008_01_var.vcf.gz
	php ../link_test_data.php 24067LRa008_01_cnvs_clincnv.tsv
	php ../link_test_data.php 24067LRa026_01.cram
	php ../link_test_data.php 24067LRa026_01.cram.crai
	php ../link_test_data.php 24067LRa026_01_var.vcf.gz
	php ../link_test_data.php 24067LRa026_01_cnvs_clincnv.tsv

call:
	php ../../src/Pipelines/multisample_longread.php -bams 24067LRa008_01.cram 24067LRa026_01.cram -status affected control -prefix test_multi -out_folder . -system system.ini -steps vc,cn,sv,an

diff:
	php ../compare_variants.php -in1 test_multi.GSvar.out -in2 test_multi.GSvar -add 24067LRa008_01,24067LRa026_01
	php ../compare_tsv.php -in1 test_multi_cnvs_clincnv.tsv.out -in2 test_multi_cnvs_clincnv.tsv -skip_cols ngsd_pathogenic_cnvs,genes,gene_info
	php ../compare_tsv.php -in1 test_multi_var_structural_variants.bedpe.out -in2 test_multi_var_structural_variants.bedpe -skip_comments_matching command=,fileDate= -skip_cols GENES,GENES_BREAKPOINTS,GENE_INFO,NGSD_HOM,NGSD_HET,NGSD_AF,NGSD_group,NGSD_SV_BREAKPOINT_DENSITY,NGSD_SV_BREAKPOINT_DENSITY_SYS,NGSD_PATHOGENIC_SVS
