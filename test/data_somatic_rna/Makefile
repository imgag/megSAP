# tumor only not tested, as it is equivalent to a single analyze_rna call
# mapping step disabled (test uses precomupted bams, counts, fusions)

X := I17R018Ra01_01
Y := I17R018Rb01_01
out_folder := project_folder

all: clean link analyze diff

clean:
	rm -rf $(out_folder)

link:
	mkdir -p $(out_folder)/Sample_$(X)
	php ../link_test_data.php $(X).bam $(out_folder)/Sample_$(X)/$(X).bam
	php ../link_test_data.php $(X).bam.bai $(out_folder)/Sample_$(X)/$(X).bam.bai
	cp -l $(X)_counts.tsv $(X)_var_fusions.tsv $(out_folder)/Sample_$(X)/
	
	mkdir -p $(out_folder)/Sample_$(Y)
	php ../link_test_data.php $(Y).bam $(out_folder)/Sample_$(Y)/$(Y).bam
	php ../link_test_data.php $(Y).bam.bai $(out_folder)/Sample_$(Y)/$(Y).bam.bai
	cp -l $(Y)_counts.tsv $(Y)_var_fusions.tsv $(out_folder)/Sample_$(Y)/
	
	mkdir -p $(out_folder)/Somatic_$(X)-$(Y)

analyze:
	php ../../src/Pipelines/somatic_rna.php \
		-p_folder $(out_folder) \
		-t_id $(X) -n_id $(Y) \
		-t_sys system.ini -n_sys system.ini \
		-o_folder $(out_folder)/Somatic_$(X)-$(Y) \
		-steps rc,fu,vc,an \
		-debug \
		--log somatic_rna_out.log

diff:
	diff $(out_folder)/Somatic_$(X)-$(Y)/$(X)-$(Y)_counts.tsv results/$(X)-$(Y)_counts.tsv
	diff $(out_folder)/Somatic_$(X)-$(Y)/$(X)-$(Y)_var_fusions.tsv results/$(X)-$(Y)_var_fusions.tsv
	php ../compare_variants.php $(out_folder)/Somatic_$(X)-$(Y)/$(X)-$(Y).GSvar results/$(X)-$(Y).GSvar
