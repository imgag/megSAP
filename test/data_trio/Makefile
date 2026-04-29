
all:
	make clear link call diff reanno
	
clear:
	rm -rf *.cram *.crai *.log *.GSvar *.vcf.gz *.vcf.gz.tbi *.qcML *_clincnv.tsv  trio_upd.tsv *.bedpe manta_evid/ reanno

link:
	php ../link_test_data.php GS140127_01.cram
	php ../link_test_data.php GS140127_01.cram.crai
	php ../link_test_data.php GS140127_01_var.vcf.gz
	php ../link_test_data.php GS140127_01_cnvs_clincnv.tsv
	php ../link_test_data.php GS140549_01.cram
	php ../link_test_data.php GS140549_01.cram.crai
	php ../link_test_data.php GS140549_01_var.vcf.gz
	php ../link_test_data.php GS140549_01_cnvs_clincnv.tsv
	php ../link_test_data.php GS140550_01.cram
	php ../link_test_data.php GS140550_01.cram.crai
	php ../link_test_data.php GS140550_01_var.vcf.gz
	php ../link_test_data.php GS140550_01_cnvs_clincnv.tsv

call:
	php ../../src/Pipelines/trio.php -c GS140127_01.cram -m GS140549_01.cram -f GS140550_01.cram -system system.ini -out_folder . -steps vc,cn,sv -no_check

diff:
	#variants
	php ../compare_variants.php -in1 trio.GSvar.out -in2 trio.GSvar -add GS140127_01,GS140549_01,GS140550_01
	#CNVs
	php ../compare_tsv.php -in1 trio_cnvs_clincnv.tsv.out -in2 trio_cnvs_clincnv.tsv -skip_cols ngsd_pathogenic_cnvs,genes,gene_info
	#UPDs
	diff -b trio_upd.tsv.out trio_upd.tsv
	#SVs
	php ../compare_tsv.php -in1 trio_var_structural_variants.bedpe.out -in2 trio_var_structural_variants.bedpe -skip_comments_matching cmdline=,fileDate= -skip_cols GENES,GENES_BREAKPOINTS,GENE_INFO,NGSD_HOM,NGSD_HET,NGSD_AF,NGSD_group,NGSD_SV_BREAKPOINT_DENSITY,NGSD_SV_BREAKPOINT_DENSITY_SYS,NGSD_PATHOGENIC_SVS

reanno:
	mkdir -p reanno
	# copy files
	cp trio_var.vcf.gz reanno/
	cp trio_cnvs_clincnv.tsv reanno/
	cp trio_var_structural_variants.bedpe reanno/
	# run reannotation
	php ../../src/Pipelines/trio.php -c GS140127_01.cram -m GS140549_01.cram -f GS140550_01.cram -system system.ini -out_folder reanno -steps vc,cn,sv -annotation_only
	# diff
	diff trio_cnvs_clincnv.tsv reanno/trio_cnvs_clincnv.tsv
	diff trio.GSvar reanno/trio.GSvar
	diff trio_var_structural_variants.bedpe reanno/trio_var_structural_variants.bedpe

######################### AUXILARY TARGETS #########################

update_cnv_data:
	php ../../src/Pipelines/analyze.php -steps cn -name GS140127_01 -folder /mnt/storage2/projects/diagnostic/KM/Sample_GS140127_01/
	php ../../src/Pipelines/analyze.php -steps cn -name GS140549_01 -folder /mnt/storage2/projects/diagnostic/KM/Sample_GS140549_01/
	php ../../src/Pipelines/analyze.php -steps cn -name GS140550_01 -folder /mnt/storage2/projects/diagnostic/KM/Sample_GS140550_01/
	cp /mnt/storage2/projects/diagnostic/KM/Sample_GS140127_01/GS140127_01_cnvs_clincnv.tsv .
	cp /mnt/storage2/projects/diagnostic/KM/Sample_GS140549_01/GS140549_01_cnvs_clincnv.tsv .
	cp /mnt/storage2/projects/diagnostic/KM/Sample_GS140550_01/GS140550_01_cnvs_clincnv.tsv .
	cp /mnt/storage2/projects/diagnostic/KM/Sample_GS140127_01/GS140127_01_cnvs_clincnv.tsv ../data_multi/
	cp /mnt/storage2/projects/diagnostic/KM/Sample_GS140549_01/GS140549_01_cnvs_clincnv.tsv ../data_multi/
	cp /mnt/storage2/projects/diagnostic/KM/Sample_GS140550_01/GS140550_01_cnvs_clincnv.tsv ../data_multi/
