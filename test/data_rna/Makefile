
all: clean analyze_se diff_se analyze_pe diff_pe

clean:
	rm -rf *.log *.txt single_end paired_end
	
analyze_se:
	mkdir -p single_end
	php ../../src/Pipelines/analyze_rna.php -in_for analyze_rna_in5.fastq.gz analyze_rna_in6.fastq.gz -out_folder single_end -out_name E16R004a01_01 -keepUnmapped -noIndelRealign -gtfFile /mnt/share/data/genomes/Ecoli/W3110/Annotation/UCSC_refSeq_ncRNA_sRNA_exonOnly_uniqSource.gtf -reference /mnt/share/data/genomes/Ecoli/W3110/STAR -steps ma,rc -system analyze_rna_ini2.ini --log analyze_rna_out2.log 

diff_se:
	diff analyze_rna_out2.tsv single_end/E16R004a01_01_fpkm.tsv

analyze_pe:
	mkdir -p paired_end
	php ../../src/Pipelines/analyze_rna.php -in_for analyze_rna_in1.fastq.gz analyze_rna_in2.fastq.gz -in_rev analyze_rna_in3.fastq.gz analyze_rna_in4.fastq.gz -out_folder paired_end -out_name I16R003a01_01 -steps ma,rc,fu,an -noIndelRealign -system analyze_rna_ini1.ini --log analyze_rna_out1.log
	
diff_pe:
	diff analyze_rna_out1.tsv paired_end/I16R003a01_01_fpkm.tsv

#################################### DEBUG TARGETS ####################################
	
fast: clean analyze_pe diff_pe
