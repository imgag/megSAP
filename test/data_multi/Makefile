
all:
	make clear link vc cn diff ra
	
clear:
	rm -rf *.bam *.bai *.log *.GSvar *.vcf.gz *.vcf.gz.tbi *.qcML trio_cnvs_clincnv.tsv

link:
	php ../link_test_data.php GS140127_01.bam
	php ../link_test_data.php GS140127_01.bam.bai
	php ../link_test_data.php GS140549_01.bam
	php ../link_test_data.php GS140549_01.bam.bai
	php ../link_test_data.php GS140550_01.bam
	php ../link_test_data.php GS140550_01.bam.bai

vc:
	rm -rf multi_log1_vc.log
	php ../../src/Pipelines/multisample.php -bams GS140127_01.bam GS140549_01.bam GS140550_01.bam -status affected control control -system system.ini -out_folder . --log multi_log1_vc.log -steps vc

cn:
	rm -rf multi_log3_cn.log
	php ../../src/Pipelines/multisample.php -bams GS140127_01.bam GS140549_01.bam GS140550_01.bam -status affected control control -system system.ini -out_folder . --log multi_log3_cn.log -steps cn

diff:
	#variants
	php ../compare_variants.php -in1 multi.GSvar.out -in2 multi.GSvar
	#CNVs
	diff multi_cnvs_clincnv.tsv multi_cnvs_clincnv.tsv.out

ra:
	mkdir -p reanno
	# copy files
	cp all.vcf.gz reanno/
	cp all_mito.vcf.gz reanno/ 
	cp multi_cnvs_clincnv.tsv reanno/
	# run reannotation
	php ../../src/Pipelines/multisample.php -bams GS140127_01.bam GS140549_01.bam GS140550_01.bam -status affected control control -system system.ini -out_folder reanno --log multi_reanno.log -steps vc,cn -annotation_only
	# diff
	diff multi_cnvs_clincnv.tsv reanno/multi_cnvs_clincnv.tsv
	diff multi.GSvar reanno/multi.GSvar

