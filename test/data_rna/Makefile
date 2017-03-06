
all: clean analyze_se diff_se analyze_pe diff_pe

clean:
	rm -rf *.log single_end paired_end
	
analyze_se:
	mkdir -p single_end
	php ../../src/Pipelines/analyze_rna.php -in_for analyze_rna_in1.fastq.gz analyze_rna_in2.fastq.gz -out_folder single_end -out_name I16R003a01_01 -steps ma,rc,an -system system.ini --log analyze_rna_out_se.log

diff_se:
	diff analyze_rna_out_se.tsv single_end/I16R003a01_01_counts_fpkm.tsv

analyze_pe:
	mkdir -p paired_end
	php ../../src/Pipelines/analyze_rna.php -in_for analyze_rna_in1.fastq.gz analyze_rna_in2.fastq.gz -in_rev analyze_rna_in3.fastq.gz analyze_rna_in4.fastq.gz -out_folder paired_end -out_name I16R003a01_01 -steps ma,rc,an,fu -system system.ini --log analyze_rna_out_pe.log
	
diff_pe:
	diff analyze_rna_out_pe.tsv paired_end/I16R003a01_01_counts_fpkm.tsv

#################################### DEBUG TARGETS ####################################
	
fast: clean analyze_pe diff_pe

fusion:
	mkdir -p paired_end
	php ../../src/Pipelines/analyze_rna.php -in_for analyze_rna_in1.fastq.gz analyze_rna_in2.fastq.gz -in_rev analyze_rna_in3.fastq.gz analyze_rna_in4.fastq.gz -out_folder paired_end -out_name I16R003a01_01 -steps fu -system system.ini --log analyze_rna_out_pe.log
	
