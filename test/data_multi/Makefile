
all:
	make clear link call diff reanno
	
clear:
	rm -rf *.cram *.crai *.log *.GSvar *.vcf.gz *.vcf.gz.tbi *.qcML multi_cnvs_clincnv.tsv *.bedpe manta_evid/ reanno 

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
	php ../../src/Pipelines/multisample.php -bams GS140127_01.cram GS140549_01.cram GS140550_01.cram -status affected control control -system system.ini -out_folder . -steps vc,cn,sv
	# remove gene annotation
	cut -f1-24 multi_var_structural_variants.bedpe > multi_var_structural_variants_cut.bedpe
	# create file from OMIM annotation and CNV overlap
	cat multi_var_structural_variants.bedpe | egrep -v "^##"| rev | cut -f 1,2| rev > multi_var_structural_variants_cut2.bedpe

diff:
	#variants
	php ../compare_variants.php -in1 multi.GSvar.out -in2 multi.GSvar -add GS140127_01,GS140549_01,GS140550_01
	#CNVs
	php ../compare_tsv.php -in1 multi_cnvs_clincnv.tsv.out -in2 multi_cnvs_clincnv.tsv -skip_cols ngsd_pathogenic_cnvs,genes,gene_info
	#SVs
	diff multi_var_structural_variants_cut.bedpe multi_var_structural_variants_cut.bedpe.out --ignore-matching-lines="^##"
	diff multi_var_structural_variants_cut2.bedpe multi_var_structural_variants_cut2.bedpe.out --ignore-matching-lines="^##"

reanno:
	mkdir -p reanno
	# copy files
	cp all.vcf.gz reanno/
	cp all_mito.vcf.gz reanno/ 
	cp multi_cnvs_clincnv.tsv reanno/
	cp multi_var_structural_variants.bedpe reanno/
	# run reannotation
	php ../../src/Pipelines/multisample.php -bams GS140127_01.cram GS140549_01.cram GS140550_01.cram -status affected control control -system system.ini -out_folder reanno -steps vc,cn,sv -annotation_only
	# diff
	diff multi_cnvs_clincnv.tsv reanno/multi_cnvs_clincnv.tsv
	diff multi.GSvar reanno/multi.GSvar
	diff multi_var_structural_variants.bedpe reanno/multi_var_structural_variants.bedpe --ignore-matching-lines="^##"
