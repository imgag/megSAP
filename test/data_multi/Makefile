
all:
	make clear link call diff reanno
	
clear:
	rm -rf *.cram *.crai *.log *.GSvar *.vcf.gz *.vcf.gz.tbi *.qcML *_cnvs_clincnv.tsv *.bedpe manta_evid/ reanno 

link:
	php ../link_test_data.php GS140549_01.cram
	php ../link_test_data.php GS140549_01.cram.crai
	php ../link_test_data.php GS140549_01_var.vcf.gz
	php ../link_test_data.php GS140549_01_cnvs_clincnv.tsv
	php ../link_test_data.php GS140550_01.cram
	php ../link_test_data.php GS140550_01.cram.crai
	php ../link_test_data.php GS140550_01_var.vcf.gz
	php ../link_test_data.php GS140550_01_cnvs_clincnv.tsv

call:
	php ../../src/Pipelines/multisample.php -bams GS140549_01.cram GS140550_01.cram -status affected affected -system system.ini -out_folder . -steps vc,cn,sv

diff:
	#variants
	php ../compare_variants.php -in1 multi.GSvar.out -in2 multi.GSvar -add GS140549_01,GS140550_01
	#CNVs
	php ../compare_tsv.php -in1 multi_cnvs_clincnv.tsv.out -in2 multi_cnvs_clincnv.tsv -skip_cols ngsd_pathogenic_cnvs,genes,gene_info
	#SVs
	php ../compare_tsv.php -in1 multi_var_structural_variants.bedpe.out -in2 multi_var_structural_variants.bedpe -skip_comments_matching cmdline=,fileDate=,PIPELINE= -skip_cols GENES,GENES_BREAKPOINTS,GENE_INFO,NGSD_HOM,NGSD_HET,NGSD_AF,NGSD_group,NGSD_SV_BREAKPOINT_DENSITY,NGSD_SV_BREAKPOINT_DENSITY_SYS,NGSD_PATHOGENIC_SVS

reanno:
	mkdir -p reanno
	# copy files
	cp multi_var.vcf.gz reanno/
	cp multi_cnvs_clincnv.tsv reanno/
	cp multi_var_structural_variants.bedpe reanno/
	# run reannotation
	php ../../src/Pipelines/multisample.php -bams GS140549_01.cram GS140550_01.cram -status affected affected -system system.ini -out_folder reanno -steps vc,cn,sv -annotation_only
	# diff
	diff multi_cnvs_clincnv.tsv reanno/multi_cnvs_clincnv.tsv
	diff multi.GSvar reanno/multi.GSvar
	diff multi_var_structural_variants.bedpe reanno/multi_var_structural_variants.bedpe --ignore-matching-lines="^##"
