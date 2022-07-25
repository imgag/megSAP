all: clean prepare analyze_pe diff_pe

clean:
	rm -rf *.log paired_end

prepare:
	mkdir -p paired_end
	ln -s ../analyze_rna_in1.fastq.gz paired_end/PE_L001_R1_001.fastq.gz
	ln -s ../analyze_rna_in2.fastq.gz paired_end/PE_L002_R1_001.fastq.gz
	ln -s ../analyze_rna_in3.fastq.gz paired_end/PE_L001_R2_001.fastq.gz
	ln -s ../analyze_rna_in4.fastq.gz paired_end/PE_L002_R2_001.fastq.gz

analyze_pe:
	php ../../src/Pipelines/analyze_rna.php \
			-folder paired_end -name PE \
			-steps ma,rc,an,fu -system system.ini \
			-library_type unstranded \
			-skip_dna_reannotation \
			--log analyze_rna_out_pe.log

diff_pe:
	diff analyze_rna_out_pe.tsv paired_end/PE_counts.tsv
	diff analyze_rna_out_pe_fusion.tsv paired_end/PE_fusions_arriba.tsv
	# remove lines which are only present if NGSD can be accessed
	sed -i '/aberrant spliced gene count/d' paired_end/PE_stats_RNA.qcML
	sed -i '/covered gene count/d' paired_end/PE_stats_RNA.qcML
	php ../compare_qcml_xml.php analyze_rna_pe_stats_RNA.qcML paired_end/PE_stats_RNA.qcML
