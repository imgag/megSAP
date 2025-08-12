
all: clean link map vc cn sv re diff reanno clean2 link2 dragen_analysis diff2

clean:
	date
	rm -rf *.fastq.gz *.log *.txt *.GSvar *.vcf.gz *.vcf.gz.tbi *.cram *.crai *.qcML *.bed *_cnvs_clincnv.seg *_cnvs_clincnv.tsv *.bedpe manta_evid repeat_expansions *_bafs.igv *_circos.png *_mosaic_cnvs.tsv *_repeats_expansionhunter.json *_repeats_expansionhunter.vcf *_rohs.tsv *_ancestry.tsv dragen_variant_calls

link:
	date
	php ../link_test_data.php NA12878_17_GTAGAGGATATCCTCT_L001_R1_001.fastq.gz
	php ../link_test_data.php NA12878_17_GTAGAGGATATCCTCT_L001_R2_001.fastq.gz

map:
	date
	php ../../src/Pipelines/analyze.php -folder . -name NA12878_17 -system system.ini -steps ma -threads 1

vc:
	date
	php ../../src/Pipelines/analyze.php -folder . -name NA12878_17 -system system.ini -steps vc

cn:
	date
	php ../../src/Pipelines/analyze.php -folder . -name NA12878_17 -system system.ini -steps cn

sv:
	date
	php ../../src/Pipelines/analyze.php -folder . -name NA12878_17 -system system.ini -steps sv
	# remove gene annotation
	cut -f1-22 NA12878_17_var_structural_variants.bedpe > NA12878_17_var_structural_variants_cut.bedpe
	# create file from OMIM annotation and CNV overlap
	cat NA12878_17_var_structural_variants.bedpe | egrep -v "^##"| rev | cut -f 1,2| rev > NA12878_17_var_structural_variants_cut2.bedpe

re:
	date
	php ../../src/Pipelines/analyze.php -folder . -name NA12878_17 -system system.ini -steps re

diff:
	date
	php ../compare_variants.php -in1 NA12878_17.GSvar.out -in2 NA12878_17.GSvar -add NA12878_17
	php ../compare_qcml_xml.php NA12878_17_stats_fastq.qcML.out NA12878_17_stats_fastq.qcML
	php ../compare_qcml_xml.php NA12878_17_stats_map.qcML.out NA12878_17_stats_map.qcML 0.01
	php ../compare_qcml_xml.php NA12878_17_stats_vc.qcML.out NA12878_17_stats_vc.qcML
	php ../compare_qcml_xml.php NA12878_17_stats_other.qcML.out NA12878_17_stats_other.qcML
	diff NA12878_17_cnvs_clincnv.tsv.out NA12878_17_cnvs_clincnv.tsv --ignore-matching-lines="#"
	php ../compare_seg.php NA12878_17_cnvs_clincnv.seg.out NA12878_17_cnvs_clincnv.seg variance=0.03,value=0.05
	php ../compare_bed.php NA12878_17_nxLRRK2_SNCA_lowcov.bed.out NA12878_17_nxLRRK2_SNCA_lowcov.bed
	diff NA12878_17_var_structural_variants.bedpe.out NA12878_17_var_structural_variants_cut.bedpe --ignore-matching-lines="^##"
	diff NA12878_17_var_structural_variants_last_cols.bedpe.out NA12878_17_var_structural_variants_cut2.bedpe --ignore-matching-lines="^##"
	diff NA12878_17_repeats_expansionhunter.vcf NA12878_17_repeats_expansionhunter.vcf.out  --ignore-matching-lines="^##"
	test -f NA12878_17_circos.png

reanno:
	date
	cp NA12878_17.GSvar NA12878_17_old.GSvar
	cp NA12878_17_cnvs_clincnv.tsv NA12878_17_old_cnvs_clincnv.tsv
	cp NA12878_17_var_structural_variants.bedpe NA12878_17_old_var_structural_variants.bedpe
	php ../../src/Pipelines/analyze.php -folder . -name NA12878_17 -system system.ini -steps vc,cn,sv -annotation_only
	diff NA12878_17.GSvar NA12878_17_old.GSvar
	diff NA12878_17_cnvs_clincnv.tsv NA12878_17_old_cnvs_clincnv.tsv
	diff NA12878_17_var_structural_variants.bedpe NA12878_17_old_var_structural_variants.bedpe

clean2:
	date
	rm -rf *.vcf.gz *.vcf.gz.tbi *.cram *.crai *.qcML *.bed *_cnvs_clincnv.seg *_cnvs_clincnv.tsv *.bedpe manta_evid repeat_expansions *_bafs.igv *_circos.png *_mosaic_cnvs.tsv *_repeats_expansionhunter.json *_repeats_expansionhunter.vcf *_rohs.tsv *_ancestry.tsv dragen_variant_calls

link2: 
	date
	php ../link_test_data.php NA12878_17.cram . -copy
	php ../link_test_data.php NA12878_17.cram.crai . -copy
	mkdir -p dragen_variant_calls
	php ../link_test_data.php NA12878_17_dragen.vcf.gz dragen_variant_calls -copy
	php ../link_test_data.php NA12878_17_dragen.vcf.gz.tbi dragen_variant_calls -copy
	php ../link_test_data.php NA12878_17_dragen_svs.vcf.gz dragen_variant_calls -copy
	php ../link_test_data.php NA12878_17_dragen_svs.vcf.gz.tbi dragen_variant_calls -copy
	php ../link_test_data.php NA12878_17_dragen.gvcf.gz dragen_variant_calls -copy
	php ../link_test_data.php NA12878_17_dragen.gvcf.gz.tbi dragen_variant_calls -copy

dragen_analysis:
	date
	php ../../src/Pipelines/analyze.php -folder . -name NA12878_17 -system system.ini -steps vc,cn,sv,re -threads 5
	# remove gene annotation
	cut -f1-22 NA12878_17_var_structural_variants.bedpe > NA12878_17_var_structural_variants_cut.bedpe
	# create file from OMIM annotation and CNV overlap
	cat NA12878_17_var_structural_variants.bedpe | egrep -v "^##"| rev | cut -f 1,2| rev > NA12878_17_var_structural_variants_cut2.bedpe

diff2:
	date
	php ../compare_variants.php -in1 NA12878_17.GSvar.out2 -in2 NA12878_17.GSvar -add NA12878_17
	php ../compare_qcml_xml.php NA12878_17_stats_fastq.qcML.out NA12878_17_stats_fastq.qcML
	php ../compare_qcml_xml.php NA12878_17_stats_map.qcML.out2 NA12878_17_stats_map.qcML 0.01
	php ../compare_qcml_xml.php NA12878_17_stats_vc.qcML.out2 NA12878_17_stats_vc.qcML
	php ../compare_qcml_xml.php NA12878_17_stats_other.qcML.out2 NA12878_17_stats_other.qcML
	diff NA12878_17_cnvs_clincnv.tsv.out NA12878_17_cnvs_clincnv.tsv --ignore-matching-lines="#"
	php ../compare_seg.php NA12878_17_cnvs_clincnv.seg.out NA12878_17_cnvs_clincnv.seg variance=0.03,value=0.05
	php ../compare_bed.php NA12878_17_nxLRRK2_SNCA_lowcov.bed.out2 NA12878_17_nxLRRK2_SNCA_lowcov.bed
	diff NA12878_17_var_structural_variants.bedpe.out2 NA12878_17_var_structural_variants_cut.bedpe --ignore-matching-lines="^##"
	diff NA12878_17_var_structural_variants_last_cols.bedpe.out2 NA12878_17_var_structural_variants_cut2.bedpe --ignore-matching-lines="^##"
	diff NA12878_17_repeats_expansionhunter.vcf NA12878_17_repeats_expansionhunter.vcf.out  --ignore-matching-lines="^##"
	test -f NA12878_17_circos.png
	date
