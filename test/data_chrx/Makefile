
all: clean link map vc cn sv re diff reanno

clean:
	rm -rf *.fastq.gz *.log *.txt *.GSvar *.vcf.gz *.vcf.gz.tbi *.cram *.crai *.qcML *.bed *_cnvs_clincnv.seg *_cnvs_clincnv.tsv *.bedpe manta_evid repeat_expansions *_bafs.igv *_circos.png *_mosaic_cnvs.tsv *_repeats_expansionhunter.json *_repeats_expansionhunter.vcf *_rohs.tsv *_ancestry.tsv

link:
	php ../link_test_data.php NA12878_17_GTAGAGGATATCCTCT_L001_R1_001.fastq.gz
	php ../link_test_data.php NA12878_17_GTAGAGGATATCCTCT_L001_R2_001.fastq.gz

map:
	php ../../src/Pipelines/analyze.php -folder . -name NA12878_17 -system system.ini -steps ma -threads 1

vc:
	php ../../src/Pipelines/analyze.php -folder . -name NA12878_17 -system system.ini -steps vc

cn:
	php ../../src/Pipelines/analyze.php -folder . -name NA12878_17 -system system.ini -steps cn

sv:
	php ../../src/Pipelines/analyze.php -folder . -name NA12878_17 -system system.ini -steps sv
	# remove gene annotation
	cut -f1-22 NA12878_17_manta_var_structural.bedpe > NA12878_17_manta_var_structural_cut.bedpe
	# create file from OMIM annotation and CNV overlap
	cat NA12878_17_manta_var_structural.bedpe | egrep -v "^##"| rev | cut -f 1,2| rev > NA12878_17_manta_var_structural_cut2.bedpe

re:
	php ../../src/Pipelines/analyze.php -folder . -name NA12878_17 -system system.ini -steps re

diff:
	php ../compare_variants.php -in1 NA12878_17.GSvar.out -in2 NA12878_17.GSvar
	php ../compare_qcml_xml.php NA12878_17_stats_fastq.qcML.out NA12878_17_stats_fastq.qcML
	php ../compare_qcml_xml.php NA12878_17_stats_map.qcML.out NA12878_17_stats_map.qcML 0.01
	php ../compare_qcml_xml.php NA12878_17_stats_vc.qcML.out NA12878_17_stats_vc.qcML
	php ../compare_qcml_xml.php NA12878_17_stats_other.qcML.out NA12878_17_stats_other.qcML
	diff NA12878_17_cnvs_clincnv.tsv.out NA12878_17_cnvs_clincnv.tsv --ignore-matching-lines="#"
	numdiff NA12878_17_cnvs_clincnv.seg.out NA12878_17_cnvs_clincnv.seg --absolute-tolerance 0.04
	php ../compare_bed.php NA12878_17_nxLRRK2_SNCA_lowcov.bed.out NA12878_17_nxLRRK2_SNCA_lowcov.bed
	diff NA12878_17_manta_var_structural.bedpe.out NA12878_17_manta_var_structural_cut.bedpe --ignore-matching-lines="^##"
	diff NA12878_17_manta_var_structural_last_cols.bedpe.out NA12878_17_manta_var_structural_cut2.bedpe --ignore-matching-lines="^##"
	diff NA12878_17_repeats_expansionhunter.vcf NA12878_17_repeats_expansionhunter.vcf.out  --ignore-matching-lines="^##"
	test -f NA12878_17_circos.png

reanno:
	cp NA12878_17.GSvar NA12878_17_old.GSvar
	cp NA12878_17_cnvs_clincnv.tsv NA12878_17_old_cnvs_clincnv.tsv
	cp NA12878_17_manta_var_structural.bedpe NA12878_17_old_manta_var_structural.bedpe
	php ../../src/Pipelines/analyze.php -folder . -name NA12878_17 -system system.ini -steps vc,cn,sv -annotation_only
	diff NA12878_17.GSvar NA12878_17_old.GSvar
	diff NA12878_17_cnvs_clincnv.tsv NA12878_17_old_cnvs_clincnv.tsv
	diff NA12878_17_manta_var_structural.bedpe NA12878_17_old_manta_var_structural.bedpe
