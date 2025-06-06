
all:
	make clear link call diff reanno
	
clear:
	rm -rf *.cram *.crai *.log *.GSvar *.vcf.gz *.vcf.gz.tbi *.qcML trio_cnvs_clincnv.tsv  trio_upd.tsv *.bedpe manta_evid/ reanno

link:
	php ../link_test_data.php GS140127_01.cram
	php ../link_test_data.php GS140127_01.cram.crai
	php ../link_test_data.php GS140127_01_cnvs_clincnv.tsv
	php ../link_test_data.php GS140549_01.cram
	php ../link_test_data.php GS140549_01.cram.crai
	php ../link_test_data.php GS140549_01_cnvs_clincnv.tsv
	php ../link_test_data.php GS140550_01.cram
	php ../link_test_data.php GS140550_01.cram.crai
	php ../link_test_data.php GS140550_01_cnvs_clincnv.tsv

call:
	php ../../src/Pipelines/trio.php -c GS140127_01.cram -m GS140549_01.cram -f GS140550_01.cram -system system.ini -out_folder . -steps vc,cn,sv -no_check
	# remove gene annotation
	cut -f1-24 trio_var_structural_variants.bedpe > trio_var_structural_variants_cut.bedpe
	# create file from OMIM annotation and CNV overlap
	cat trio_var_structural_variants.bedpe | egrep -v "^##"| rev | cut -f 1,2| rev > trio_var_structural_variants_cut2.bedpe	

diff:
	#variants
	php ../compare_variants.php -in1 trio.GSvar.out -in2 trio.GSvar
	#CNVs
	php ../compare_tsv.php -in1 trio_cnvs_clincnv.tsv.out -in2 trio_cnvs_clincnv.tsv -skip_cols ngsd_pathogenic_cnvs,genes,gene_info
	#UPDs
	diff -b trio_upd.tsv.out trio_upd.tsv
	#SVs
	diff trio_var_structural_variants_cut.bedpe trio_var_structural_variants_cut.bedpe.out --ignore-matching-lines="^##"
	diff trio_var_structural_variants_cut2.bedpe trio_var_structural_variants_cut2.bedpe.out --ignore-matching-lines="^##"

reanno:
	mkdir -p reanno
	# copy files
	cp all.vcf.gz reanno/
	cp all_mito.vcf.gz reanno/
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
