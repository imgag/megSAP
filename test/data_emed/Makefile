#pipeline test emed
all: clean link emed diff finished

clean:
	rm -f *.out.vcf.gz *.out.GSvar *.out.tsv
	rm -rf Sample_GS160540_01
	rm -rf Sample_GS160539_01
	rm -rf Sample_I16R003a06_01
	rm -rf Sample_I16R003a05_01
	rm -rf results_emed1
	rm -rf results_emed2

link:
	#sample folders
	mkdir -p Sample_GS160540_01
	mkdir -p Sample_GS160539_01
	mkdir -p Sample_I16R003a06_01
	mkdir -p Sample_I16R003a05_01
	#result folders
	mkdir -p results_emed1
	mkdir -p results_emed2
	#raw data
	php ../link_test_data.php GS160540_01_*.fastq.gz Sample_GS160540_01
	php ../link_test_data.php GS160539_01*.fastq.gz Sample_GS160539_01
	php ../link_test_data.php I16R003a06_01*.fastq.gz Sample_I16R003a06_01
	php ../link_test_data.php I16R003a05_01*.fastq.gz Sample_I16R003a05_01
	#result files
	php ../link_test_data.php GS160540_01-GS160539_01*emed1.vcf.gz
	php ../link_test_data.php GS160540_01-GS160539_01*emed2.vcf.gz
	php ../link_test_data.php GS160540_01-GS160539_01*emed1.tsv
	php ../link_test_data.php GS160540_01-GS160539_01*emed2.tsv
	php ../link_test_data.php GS160540_01-GS160539_01*emed1.GSvar
	php ../link_test_data.php GS160540_01-GS160539_01*emed2.GSvar
	php ../link_test_data.php GS160539_01*emed1.vcf.gz
	php ../link_test_data.php GS160539_01*emed2.vcf.gz
	php ../link_test_data.php I16R003a06_01-I16R003a05_01*emed*.tsv

emed:
	#sample pair, rna pair, start with mapping
	php ../../src/Pipelines/somatic_emed.php \
		-p_folder . -t_dna_id GS160540_01 -n_dna_id GS160539_01 -o_folder results_emed1 \
		-t_rna_id I16R003a06_01 -t_rna_fo Sample_I16R003a06_01 -n_rna_id I16R003a05_01 -n_rna_fo Sample_I16R003a05_01 \
		-t_dna_sys system_dna.ini -n_dna_sys system_dna.ini -t_rna_sys system_rna.ini -n_rna_sys system_rna.ini \
		-steps ma,fu,vc,an,co --log results_emed1/somatic.emed1.log
	#somatic pair, only tumor rna, start with annotation
	php ../../src/Pipelines/somatic_emed.php \
		-p_folder . -t_dna_id GS160540_01 -n_dna_id GS160539_01 -o_folder results_emed2 \
		-t_rna_id I16R003a06_01 -t_rna_fo Sample_I16R003a06_01 \
		-t_dna_sys system_dna.ini -n_dna_sys system_dna.ini -t_rna_sys system_rna.ini -n_rna_sys system_rna.ini \
		-steps vc,an,co --log results_emed2/somatic.emed2.log

diff:
	#emed1
	php ../compare_variants.php GS160540_01-GS160539_01.emed1.GSvar results_emed1/GS160540_01-GS160539_01.GSvar
	php ../compare_variants.php GS160540_01-GS160539_01_var_annotated.emed1.vcf.gz results_emed1/GS160540_01-GS160539_01_var_annotated.vcf.gz
	php ../compare_variants.php GS160539_01_adme_var_annotated.emed1.vcf.gz results_emed1/GS160539_01_adme_var_annotated.vcf.gz
	diff GS160540_01-GS160539_01_var_copy.emed1.tsv results_emed1/GS160540_01-GS160539_01_cnvs.tsv
	diff I16R003a06_01-I16R003a05_01_var_fusions.emed1.tsv results_emed1/I16R003a06_01-I16R003a05_01_var_fusions.tsv
	#structural variants currently not implemented / no error message by SampleDiff	php ../compare_variants.php GS160540_01-GS160539_01_var_structural.emed1.out.vcf.gz results_emed1/GS160540_01-GS160539_01_var_structural.vcf.gz
	#emed2
	php ../compare_variants.php GS160540_01-GS160539_01.emed2.GSvar results_emed2/GS160540_01-GS160539_01.GSvar
	php ../compare_variants.php GS160540_01-GS160539_01_var_annotated.emed2.vcf.gz results_emed2/GS160540_01-GS160539_01_var_annotated.vcf.gz
	php ../compare_variants.php GS160539_01_adme_var_annotated.emed2.vcf.gz results_emed2/GS160539_01_adme_var_annotated.vcf.gz
	diff GS160540_01-GS160539_01_var_copy.emed2.tsv results_emed2/GS160540_01-GS160539_01_cnvs.tsv

finished:
	@echo "Finished pipeline test."

update_fusions:
	php ../update_test_data.php results_emed1/I16R003a06_01-I16R003a05_01_var_fusions.tsv I16R003a06_01-I16R003a05_01_var_fusions.emed1.tsv

update_variants:
	php ../update_test_data.php results_emed1/GS160540_01-GS160539_01.GSvar GS160540_01-GS160539_01.emed1.GSvar
	php ../update_test_data.php results_emed1/GS160540_01-GS160539_01_var_annotated.vcf.gz GS160540_01-GS160539_01_var_annotated.emed1.vcf.gz
	php ../update_test_data.php results_emed1/GS160539_01_adme_var_annotated.vcf.gz GS160539_01_adme_var_annotated.emed1.vcf.gz
	php ../update_test_data.php results_emed2/GS160540_01-GS160539_01.GSvar GS160540_01-GS160539_01.emed2.GSvar
	php ../update_test_data.php results_emed2/GS160540_01-GS160539_01_var_annotated.vcf.gz GS160540_01-GS160539_01_var_annotated.emed2.vcf.gz
	php ../update_test_data.php results_emed2/GS160539_01_adme_var_annotated.vcf.gz GS160539_01_adme_var_annotated.emed2.vcf.gz

update_cnvs:
	php ../update_test_data.php results_emed1/GS160540_01-GS160539_01_cnvs.tsv GS160540_01-GS160539_01_var_copy.emed1.tsv
	php ../update_test_data.php results_emed2/GS160540_01-GS160539_01_cnvs.tsv GS160540_01-GS160539_01_var_copy.emed2.tsv