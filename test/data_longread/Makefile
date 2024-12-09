
all: clean link map vc cn sv re me an diff bam_analysis diff2

clean:
	rm -rf *.fastq.gz *.log *.txt *.GSvar *.gvcf.gz *.gvcf.gz.tbi *.vcf.gz *.vcf.gz.tbi *.bam *.bai *.qcML *.bed *_cnvs_clincnv.seg *_cnvs_clincnv.tsv *.bedpe clair_temp repeat_expansions *_bafs.igv *_circos.png *_mosaic_cnvs.tsv *_rohs.tsv *_ancestry.tsv *_methylation.tsv

link:
	php ../link_test_data.php 23014LRa023L2_01_chr1_longread.fastq.gz

map:
	php ../../src/Pipelines/analyze_longread.php -folder . -name 23014LRa023L2_01 -system system.ini -steps ma -threads 1
	
vc:
	php ../../src/Pipelines/analyze_longread.php -folder . -name 23014LRa023L2_01 -system system.ini -steps vc 
	
cn:
	php ../../src/Pipelines/analyze_longread.php -folder . -name 23014LRa023L2_01 -system system.ini -steps cn

sv:
	php ../../src/Pipelines/analyze_longread.php -folder . -name 23014LRa023L2_01 -system system.ini -steps sv

re:
	php ../../src/Pipelines/analyze_longread.php -folder . -name 23014LRa023L2_01 -system system.ini -steps re

me:
	# should only print waring since input data doesn't contain methylation
	php ../../src/Pipelines/analyze_longread.php -folder . -name 23014LRa023L2_01 -system system.ini -steps me 

an:
	php ../../src/Pipelines/analyze_longread.php -folder . -name 23014LRa023L2_01 -system system.ini -steps an
	# remove gene annotation
	cut -f1-22 23014LRa023L2_01_var_structural_variants.bedpe > 23014LRa023L2_01_var_structural_variants_cut.bedpe

diff:
	php ../compare_variants.php -in1 23014LRa023L2_01.GSvar.out -in2 23014LRa023L2_01.GSvar
	php ../compare_qcml_xml.php 23014LRa023L2_01_stats_fastq.qcML.out 23014LRa023L2_01_stats_fastq.qcML
	php ../compare_qcml_xml.php 23014LRa023L2_01_stats_map.qcML.out 23014LRa023L2_01_stats_map.qcML
	php ../compare_qcml_xml.php 23014LRa023L2_01_stats_vc.qcML.out 23014LRa023L2_01_stats_vc.qcML 1.01
	php ../compare_qcml_xml.php 23014LRa023L2_01_stats_other.qcML.out 23014LRa023L2_01_stats_other.qcML
	diff 23014LRa023L2_01_cnvs_clincnv.tsv.out 23014LRa023L2_01_cnvs_clincnv.tsv --ignore-matching-lines="#"
	diff 23014LRa023L2_01_var_structural_variants.bedpe.out 23014LRa023L2_01_var_structural_variants_cut.bedpe --ignore-matching-lines="^##"
	diff 23014LRa023L2_01_prs.tsv.out 23014LRa023L2_01_prs.tsv
	diff 23014LRa023L2_01_phasing_track.bed.out 23014LRa023L2_01_phasing_track.bed
	diff 23014LRa023L2_01_bafs.igv.out 23014LRa023L2_01_bafs.igv

# 2nd test with BAM input

bam_analysis: clean
	php ../link_test_data.php 23014LRa023L2_02.mod.unmapped.bam
	php ../../src/Pipelines/analyze_longread.php -folder . -name 23014LRa023L2_02 -system system.ini -steps ma,vc,cn,sv,re,me,an
	cut -f1-22 23014LRa023L2_02_var_structural_variants.bedpe > 23014LRa023L2_02_var_structural_variants_cut.bedpe

diff2:
	php ../compare_variants.php -in1 23014LRa023L2_02.GSvar.out -in2 23014LRa023L2_02.GSvar
	php ../compare_qcml_xml.php 23014LRa023L2_02_stats_fastq.qcML.out 23014LRa023L2_02_stats_fastq.qcML
	php ../compare_qcml_xml.php 23014LRa023L2_02_stats_map.qcML.out 23014LRa023L2_02_stats_map.qcML
	php ../compare_qcml_xml.php 23014LRa023L2_02_stats_vc.qcML.out 23014LRa023L2_02_stats_vc.qcML
	php ../compare_qcml_xml.php 23014LRa023L2_02_stats_other.qcML.out 23014LRa023L2_02_stats_other.qcML
	diff 23014LRa023L2_02_cnvs_clincnv.tsv.out 23014LRa023L2_02_cnvs_clincnv.tsv --ignore-matching-lines="#"
	diff 23014LRa023L2_02_var_structural_variants.bedpe.out 23014LRa023L2_02_var_structural_variants_cut.bedpe --ignore-matching-lines="^##"
	diff 23014LRa023L2_02_prs.tsv.out 23014LRa023L2_02_prs.tsv
	diff 23014LRa023L2_02_phasing_track.bed.out 23014LRa023L2_02_phasing_track.bed
	diff 23014LRa023L2_02_bafs.igv.out 23014LRa023L2_02_bafs.igv
	diff 23014LRa023L2_02_var_methylation.tsv.out 23014LRa023L2_02_var_methylation.tsv
