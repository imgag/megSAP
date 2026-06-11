
all: default_pipeline dragen_analysis

### non-DRAGEN pipeline ###
default_pipeline: clean link map call diff reanno

clean:
	rm -rf *.fastq.gz *.log *.txt *.GSvar *.vcf.gz *.vcf.gz.tbi *.cram *.crai *.bam *.bai *.qcML *.bed *_cnvs_clincnv.seg *_cnvs_clincnv.tsv *.bedpe manta_evid repeat_expansions *_bafs.igv *_circos.png *_mosaic_cnvs.tsv *_repeats_expansionhunter.json *_repeats_expansionhunter.vcf *_rohs.tsv *_ancestry.tsv dragen

link:
	php ../link_test_data.php NA12878_17_GTAGAGGATATCCTCT_L001_R1_001.fastq.gz
	php ../link_test_data.php NA12878_17_GTAGAGGATATCCTCT_L001_R2_001.fastq.gz

map:
	php ../../src/Pipelines/analyze.php -folder . -name NA12878_17 -system system.ini -steps ma -threads 1

call:
	php ../../src/Pipelines/analyze.php -folder . -name NA12878_17 -system system.ini -steps vc,cn,sv,re

diff:
	php ../compare_variants.php -in1 NA12878_17.GSvar.out -in2 NA12878_17.GSvar -add NA12878_17
	php ../compare_qcml_xml.php NA12878_17_stats_fastq.qcML.out NA12878_17_stats_fastq.qcML
	php ../compare_qcml_xml.php NA12878_17_stats_map.qcML.out NA12878_17_stats_map.qcML 0.01
	php ../compare_qcml_xml.php NA12878_17_stats_vc.qcML.out NA12878_17_stats_vc.qcML
	php ../compare_qcml_xml.php NA12878_17_stats_other.qcML.out NA12878_17_stats_other.qcML
	diff NA12878_17_cnvs_clincnv.tsv.out NA12878_17_cnvs_clincnv.tsv --ignore-matching-lines="#"
	php ../compare_seg.php NA12878_17_cnvs_clincnv.seg.out NA12878_17_cnvs_clincnv.seg variance=0.03,value=0.05
	php ../compare_bed.php NA12878_17_nxLRRK2_SNCA_lowcov.bed.out NA12878_17_nxLRRK2_SNCA_lowcov.bed
	php ../compare_tsv.php -in1 NA12878_17_var_structural_variants.bedpe.out -in2 NA12878_17_var_structural_variants.bedpe -skip_comments_matching cmdline=,fileDate= -skip_cols GENES,GENES_BREAKPOINTS,GENE_INFO,NGSD_HOM,NGSD_HET,NGSD_AF,NGSD_group,NGSD_SV_BREAKPOINT_DENSITY,NGSD_SV_BREAKPOINT_DENSITY_SYS,NGSD_PATHOGENIC_SVS
	diff NA12878_17_repeats_expansionhunter.vcf NA12878_17_repeats_expansionhunter.vcf.out  --ignore-matching-lines="^##"
	test -f NA12878_17_circos.png

reanno:
	cp NA12878_17.GSvar NA12878_17_old.GSvar
	cp NA12878_17_cnvs_clincnv.tsv NA12878_17_old_cnvs_clincnv.tsv
	cp NA12878_17_var_structural_variants.bedpe NA12878_17_old_var_structural_variants.bedpe
	php ../../src/Pipelines/analyze.php -folder . -name NA12878_17 -system system.ini -steps vc,cn,sv -annotation_only
	diff NA12878_17.GSvar NA12878_17_old.GSvar
	diff NA12878_17_cnvs_clincnv.tsv NA12878_17_old_cnvs_clincnv.tsv
	diff NA12878_17_var_structural_variants.bedpe NA12878_17_old_var_structural_variants.bedpe

### DRAGEN pipeline ###
dragen_analysis: clean2 link2 dragen_call diff2

clean2:
	rm -rf *.vcf.gz *.vcf.gz.tbi *.cram *.crai *.qcML *.bed *_cnvs_clincnv.seg *_cnvs_clincnv.tsv *.bedpe manta_evid repeat_expansions *_bafs.igv *_circos.png *_mosaic_cnvs.tsv *_repeats_expansionhunter.json *_repeats_expansionhunter.vcf *_rohs.tsv *_ancestry.tsv dragen

link2: 
	mkdir -p dragen
	php ../link_test_data.php NA12878_17.cram dragen -copy
	php ../link_test_data.php NA12878_17.cram.crai dragen -copy
	php ../link_test_data.php NA12878_17.hard-filtered.vcf.gz dragen -copy
	php ../link_test_data.php NA12878_17.hard-filtered.vcf.gz.tbi dragen -copy
	php ../link_test_data.php NA12878_17.sv.vcf.gz dragen -copy
	php ../link_test_data.php NA12878_17.sv.vcf.gz.tbi dragen -copy
	php ../link_test_data.php NA12878_17.hard-filtered.gvcf.gz dragen -copy
	php ../link_test_data.php NA12878_17.hard-filtered.gvcf.gz.tbi dragen -copy

dragen_call:
	php ../../src/Pipelines/analyze.php -folder . -name NA12878_17 -system system.ini -steps vc,cn,sv,re -threads 5

diff2:
	php ../compare_variants.php -in1 NA12878_17.GSvar.out2 -in2 NA12878_17.GSvar -add NA12878_17
	php ../compare_qcml_xml.php NA12878_17_stats_fastq.qcML.out NA12878_17_stats_fastq.qcML
	php ../compare_qcml_xml.php NA12878_17_stats_map.qcML.out2 NA12878_17_stats_map.qcML 0.01
	php ../compare_qcml_xml.php NA12878_17_stats_vc.qcML.out2 NA12878_17_stats_vc.qcML
	php ../compare_qcml_xml.php NA12878_17_stats_other.qcML.out2 NA12878_17_stats_other.qcML
	diff NA12878_17_cnvs_clincnv.tsv.out NA12878_17_cnvs_clincnv.tsv --ignore-matching-lines="#"
	php ../compare_seg.php NA12878_17_cnvs_clincnv.seg.out NA12878_17_cnvs_clincnv.seg variance=0.03,value=0.05
	php ../compare_bed.php NA12878_17_nxLRRK2_SNCA_lowcov.bed.out2 NA12878_17_nxLRRK2_SNCA_lowcov.bed
	php ../compare_tsv.php -in1 NA12878_17_var_structural_variants.bedpe.out2 -in2 NA12878_17_var_structural_variants.bedpe -skip_comments_matching cmdline=,fileDate= -skip_cols GENES,GENES_BREAKPOINTS,GENE_INFO,NGSD_HOM,NGSD_HET,NGSD_AF,NGSD_group,NGSD_SV_BREAKPOINT_DENSITY,NGSD_SV_BREAKPOINT_DENSITY_SYS,NGSD_PATHOGENIC_SVS
	diff NA12878_17_repeats_expansionhunter.vcf NA12878_17_repeats_expansionhunter.vcf.out  --ignore-matching-lines="^##"
	test -f NA12878_17_circos.png
