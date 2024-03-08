
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
	php ../link_test_data.php 21073LRa015_01.bam
	php ../link_test_data.php 21073LRa015_01.bam.bai
	php ../link_test_data.php 21073LRa015_01_cnvs_clincnv.tsv
	php ../link_test_data.php 21073LRa015_01_var.gvcf.gz
	php ../link_test_data.php 21073LRa015_01_var.gvcf.gz.tbi


vc:
	php ../../src/Pipelines/trio_longread.php -c 21073LRa015_01.bam -f 21073LRa013_01.bam -m 21073LRa014_01.bam -out_folder . -no_check -system system.ini -steps vc 
	
cn:
	php ../../src/Pipelines/trio_longread.php -c 21073LRa015_01.bam -f 21073LRa013_01.bam -m 21073LRa014_01.bam -out_folder . -no_check -system system.ini -steps cn 

sv:
	php ../../src/Pipelines/trio_longread.php -c 21073LRa015_01.bam -f 21073LRa013_01.bam -m 21073LRa014_01.bam -out_folder . -no_check -system system.ini -steps sv 

an:
	php ../../src/Pipelines/trio_longread.php -c 21073LRa015_01.bam -f 21073LRa013_01.bam -m 21073LRa014_01.bam -out_folder . -no_check -system system.ini -steps an 
	# remove gene annotation
	cut -f1-24 trio_var_structural_variants.bedpe > trio_var_structural_variants_cut.bedpe
	# create file from the constant last columns
	cat trio_var_structural_variants.bedpe | egrep -v "^##"| rev | cut -f 1-8| rev > trio_var_structural_variants_cut2.bedpe
	# remove CNV NGSD annotation
	cut -f1-18 trio_cnvs_clincnv.tsv > trio_cnvs_clincnv_cut.tsv

diff:
	php ../compare_variants.php -in1 trio.GSvar.out -in2 trio.GSvar
	diff trio_cnvs_clincnv.tsv.out trio_cnvs_clincnv_cut.tsv --ignore-matching-lines="#"
	diff trio_var_structural_variants.bedpe.out trio_var_structural_variants_cut.bedpe --ignore-matching-lines="^##"
	diff trio_var_structural_variants.bedpe.out2 trio_var_structural_variants_cut2.bedpe --ignore-matching-lines="^##"

