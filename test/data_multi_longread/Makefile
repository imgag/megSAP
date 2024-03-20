all: clean link vc cn sv an diff

clean:
	rm -rf *.fastq.gz *.log *.txt *.GSvar *.vcf.gz *.vcf.gz.tbi *.bam *.bai *.qcML *.bed *_cnvs_clincnv.seg *_cnvs_clincnv.tsv *.bedpe clair_temp repeat_expansions *_bafs.igv *_circos.png *_mosaic_cnvs.tsv *_rohs.tsv *_ancestry.tsv *.gvcf.gz *.gvcf.gz.tbi

link:
	php ../link_test_data.php 21073LRa013_01.bam
	php ../link_test_data.php 21073LRa013_01.bam.bai
	php ../link_test_data.php 21073LRa013_01_cnvs_clincnv.tsv
	php ../link_test_data.php 21073LRa013_01_var.gvcf.gz
	php ../link_test_data.php 21073LRa013_01_var.gvcf.gz.tbi
	php ../link_test_data.php 21073LRa014_01.bam
	php ../link_test_data.php 21073LRa014_01.bam.bai
	php ../link_test_data.php 21073LRa014_01_cnvs_clincnv.tsv
	php ../link_test_data.php 21073LRa014_01_var.gvcf.gz
	php ../link_test_data.php 21073LRa014_01_var.gvcf.gz.tbi


vc:
	php ../../src/Pipelines/multisample_longread.php -bams 21073LRa013_01.bam 21073LRa014_01.bam -status affected control -prefix test_multi -out_folder . -system system.ini -steps vc 
	
cn:
	php ../../src/Pipelines/multisample_longread.php -bams 21073LRa013_01.bam 21073LRa014_01.bam -status affected control -prefix test_multi -out_folder . -system system.ini -steps cn 

sv:
	php ../../src/Pipelines/multisample_longread.php -bams 21073LRa013_01.bam 21073LRa014_01.bam -status affected control -prefix test_multi -out_folder . -system system.ini -steps sv 

an:
	php ../../src/Pipelines/multisample_longread.php -bams 21073LRa013_01.bam 21073LRa014_01.bam -status affected control -prefix test_multi -out_folder . -system system.ini -steps an 
	# remove gene annotation
	cut -f1-23 test_multi_var_structural_variants.bedpe > test_multi_var_structural_variants_cut.bedpe
	# create file from OMIM annotation and CNV overlap
	cat test_multi_var_structural_variants.bedpe | egrep -v "^##"| rev | cut -f 1-7| rev > test_multi_var_structural_variants_cut2.bedpe
	# remove CNV NGSD annotation
	cut -f1-18 test_multi_cnvs_clincnv.tsv > test_multi_cnvs_clincnv_cut.tsv

diff:
	php ../compare_variants.php -in1 test_multi.GSvar.out -in2 test_multi.GSvar
	diff test_multi_cnvs_clincnv.tsv.out test_multi_cnvs_clincnv_cut.tsv --ignore-matching-lines="#"
	diff test_multi_var_structural_variants.bedpe.out test_multi_var_structural_variants_cut.bedpe --ignore-matching-lines="^##"
	diff test_multi_var_structural_variants.bedpe.out2 test_multi_var_structural_variants_cut2.bedpe --ignore-matching-lines="^##"

