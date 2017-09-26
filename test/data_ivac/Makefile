#pipeline test ivac
all: clean link ivac diff finished

clean:
	rm -f *.out.vcf.gz *.out.GSvar *.out.tsv
	rm -rf Sample_GS140850_02
	rm -rf Sample_GS140851_02
	rm -rf Sample_I14R010d01_01
	rm -rf Sample_I14R010d02_01
	rm -rf results_ivac1
	rm -rf results_ivac2

link:
	#sample folders
	mkdir -p Sample_GS140850_02
	mkdir -p Sample_GS140851_02
	mkdir -p Sample_I14R010d01_01
	mkdir -p Sample_I14R010d02_01
	#result folders
	mkdir -p results_ivac1
	mkdir -p results_ivac2
	#raw data
	php ../link_test_data.php GS140850_02_*.fastq.gz Sample_GS140850_02
	php ../link_test_data.php GS140851_02*.fastq.gz Sample_GS140851_02
	php ../link_test_data.php I14R010d01_01*.fastq.gz Sample_I14R010d01_01
	php ../link_test_data.php I14R010d02_01*.fastq.gz Sample_I14R010d02_01
	#result files
	php ../link_test_data.php GS140850_02-GS140851_02_*ivac*.out.vcf.gz .
	php ../link_test_data.php GS140850_02-GS140851_02_*ivac*.out.tsv .
	php ../link_test_data.php GS140850_02-GS140851_02_*ivac*.out.GSvar .
	php ../link_test_data.php I14R010d01_01-I14R010d02_01_*ivac*.out.tsv
	php ../link_test_data.php ivac_qcml.out.txt

ivac:
	#sample pair, rna pair, start with mapping
	php ../../src/Pipelines/somatic_ivac.php \
		-p_folder . -t_dna_id GS140850_02 -n_dna_id GS140851_02 -o_folder results_ivac1 \
		-t_rna_id I14R010d01_01 -t_rna_fo Sample_I14R010d01_01 -n_rna_id I14R010d02_01 -n_rna_fo Sample_I14R010d02_01 \
		-t_dna_sys system_dna.ini -n_dna_sys system_dna.ini -t_rna_sys system_rna.ini -n_rna_sys system_rna.ini \
		-steps ma,vc,fu,an,co,im --log results_ivac1/somatic_ivac1.log
	#somatic pair, only tumor rna, start with annotation
	php ../../src/Pipelines/somatic_ivac.php \
		-p_folder . -t_dna_id GS140850_02 -n_dna_id GS140851_02 -o_folder results_ivac2 \
		-t_rna_id I14R010d01_01 -t_rna_fo Sample_I14R010d01_01 \
		-t_dna_sys system_dna.ini -n_dna_sys system_dna.ini -t_rna_sys system_rna.ini -n_rna_sys system_rna.ini \
		-steps vc,an,co,im --log results_ivac2/somatic_ivac2.log

diff:
	#TODO: check germline vcf_annotated, add copy and structural variants for ivac2
	#ivac1
	php ../compare_variants.php GS140850_02-GS140851_02_ivac1.out.GSvar results_ivac1/GS140850_02-GS140851_02.GSvar
	php ../compare_variants.php GS140850_02-GS140851_02_var_annotated_ivac1.out.vcf.gz results_ivac1/GS140850_02-GS140851_02_var_annotated.vcf.gz
	diff GS140850_02-GS140851_02_var_copy_ivac1.out.tsv results_ivac1/GS140850_02-GS140851_02_cnvs.tsv
	diff I14R010d01_01-I14R010d02_01_var_fusions_ivac1.out.tsv results_ivac1/I14R010d01_01-I14R010d02_01_var_fusions.tsv
	#ivac2
	php ../compare_variants.php GS140850_02-GS140851_02_ivac2.out.GSvar results_ivac2/GS140850_02-GS140851_02.GSvar
	php ../compare_variants.php GS140850_02-GS140851_02_var_annotated_ivac2.out.vcf.gz results_ivac2/GS140850_02-GS140851_02_var_annotated.vcf.gz
	diff GS140850_02-GS140851_02_var_copy_ivac2.out.tsv results_ivac2/GS140850_02-GS140851_02_cnvs.tsv
	#qc
	ls */*.qcML | sort | xargs grep "<qualityParameter" > ivac_qcml.txt
	php ../compare_qcml.php ivac_qcml.out.txt ivac_qcml.txt

finished:
	@echo "Finished pipeline test."

update_variants:
	php ../update_test_data.php results_ivac1/GS140850_02-GS140851_02.GSvar GS140850_02-GS140851_02_ivac1.out.GSvar
	php ../update_test_data.php results_ivac1/GS140850_02-GS140851_02_var_annotated.vcf.gz GS140850_02-GS140851_02_var_annotated_ivac1.out.vcf.gz
	php ../update_test_data.php results_ivac2/GS140850_02-GS140851_02.GSvar GS140850_02-GS140851_02_ivac2.out.GSvar
	php ../update_test_data.php results_ivac2/GS140850_02-GS140851_02_var_annotated.vcf.gz GS140850_02-GS140851_02_var_annotated_ivac2.out.vcf.gz

update_fusions:
	php ../update_test_data.php results_ivac1/I14R010d01_01-I14R010d02_01_var_fusions.tsv I14R010d01_01-I14R010d02_01_var_fusions_ivac1.out.tsv

update_cnvs:
	php ../update_test_data.php results_ivac1/GS140850_02-GS140851_02_cnvs.tsv GS140850_02-GS140851_02_var_copy_ivac1.out.tsv
	php ../update_test_data.php results_ivac2/GS140850_02-GS140851_02_cnvs.tsv GS140850_02-GS140851_02_var_copy_ivac2.out.tsv

update_qc:
	php ../update_test_data.php ivac_qcml.txt ivac_qcml.out.txt