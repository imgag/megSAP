
all:
	make clear link vc cn sv diff ra
	
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
	php ../../src/Pipelines/multisample.php -bams GS140127_01.bam GS140549_01.bam GS140550_01.bam -status affected control control -system system.ini -out_folder . -steps vc

cn:
	php ../../src/Pipelines/multisample.php -bams GS140127_01.bam GS140549_01.bam GS140550_01.bam -status affected control control -system system.ini -out_folder . -steps cn

sv:
	php ../../src/Pipelines/multisample.php -bams GS140127_01.bam GS140549_01.bam GS140550_01.bam -status affected control control -system system.ini -out_folder . -steps sv

diff:
	#variants
	php ../compare_variants.php -in1 multi.GSvar.out -in2 multi.GSvar
	#CNVs
	diff multi_cnvs_clincnv.tsv multi_cnvs_clincnv.tsv.out
	#SVs
	diff multi_manta_var_structural.bedpe multi_manta_var_structural.bedpe.out --ignore-matching-lines="^##"

ra:
	mkdir -p reanno
	# copy files
	cp all.vcf.gz reanno/
	cp all_mito.vcf.gz reanno/ 
	cp multi_cnvs_clincnv.tsv reanno/
	cp multi_manta_var_structural.bedpe reanno/
	# run reannotation
	php ../../src/Pipelines/multisample.php -bams GS140127_01.bam GS140549_01.bam GS140550_01.bam -status affected control control -system system.ini -out_folder reanno -steps vc,cn,sv -annotation_only
	# diff
	diff multi_cnvs_clincnv.tsv reanno/multi_cnvs_clincnv.tsv
	diff multi.GSvar reanno/multi.GSvar
	diff multi_manta_var_structural.bedpe reanno/multi_manta_var_structural.bedpe
