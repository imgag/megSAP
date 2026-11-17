TIMESTAMP=$(shell date +"%Y-%m-%d_%H-%M-%S")

# only run 1 bam test to speed up tests
pipeline_test:
	make clean
	make clean_logs
	make test_bam

all: 
	make clean
	make clean_logs
	make test_fastq
	make test_bam
	make test_bam_single_steps

clean:
	rm -rf *.fastq.gz *.txt *.GSvar.gz *.gvcf.gz *.gvcf.gz.tbi *.vcf.gz *.vcf.gz.tbi *.cram *.crai *.bam *.bai *.qcML *.bed *_cnvs_clincnv.seg *_cnvs_clincnv.tsv *.bedpe clair_temp repeat_expansions *_bafs.igv *_circos.png *_mosaic_cnvs.tsv *_rohs.tsv *_ancestry.tsv  *_methylation.tsv methylartist *_modkit_track.bed.gz *_modkit_track.bed.gz.tbi *_prs.tsv *_repeats.tsv *_repeats.vcf *_var_modcall.vcf paraphase

clean_logs:
	rm -rf *.log

link:
	php ../link_test_data.php 23014LRa023L2_01_chr1_longread.fastq.gz

test_fastq:
	make clean
	make link
	make fastq_analysis
	make fastq_diff

test_bam:
	make clean
	make link
	make bam_analysis
	make bam_diff

test_bam_single_steps:
	make clean
	make link
	make bam_analysis_single_steps
	make bam_diff

# 1st test with FastQ input and performing each step separate (not performed on default)
fastq_analysis:
	php ../../src/Pipelines/analyze_longread.php -folder . -name 23014LRa023L2_01 -system system.ini -steps ma -threads 4 --log analyze_longread_fastq_ma_$(TIMESTAMP).log
	php ../../src/Pipelines/analyze_longread.php -folder . -name 23014LRa023L2_01 -system system.ini -steps vc -threads 4 --log analyze_longread_fastq_vc_$(TIMESTAMP).log
	php ../../src/Pipelines/analyze_longread.php -folder . -name 23014LRa023L2_01 -system system.ini -steps cn -threads 4 --log analyze_longread_fastq_cn_$(TIMESTAMP).log -test
	php ../../src/Pipelines/analyze_longread.php -folder . -name 23014LRa023L2_01 -system system.ini -steps sv -threads 4 --log analyze_longread_fastq_sv_$(TIMESTAMP).log
	php ../../src/Pipelines/analyze_longread.php -folder . -name 23014LRa023L2_01 -system system.ini -steps ph -threads 4 --log analyze_longread_fastq_ph_$(TIMESTAMP).log
	php ../../src/Pipelines/analyze_longread.php -folder . -name 23014LRa023L2_01 -system system.ini -steps re -threads 4 --log analyze_longread_fastq_re_$(TIMESTAMP).log
	# should only print warning since input data doesn't contain methylation
	php ../../src/Pipelines/analyze_longread.php -folder . -name 23014LRa023L2_01 -system system.ini -steps me -threads 4 --log analyze_longread_fastq_me_$(TIMESTAMP).log
	php ../../src/Pipelines/analyze_longread.php -folder . -name 23014LRa023L2_01 -system system.ini -steps pg -threads 4 --log analyze_longread_fastq_pg_$(TIMESTAMP).log
	php ../../src/Pipelines/analyze_longread.php -folder . -name 23014LRa023L2_01 -system system.ini -steps an -threads 4 --log analyze_longread_fastq_an_$(TIMESTAMP).log

fastq_diff:
	php ../compare_variants.php -in1 23014LRa023L2_01.GSvar.gz.out -in2 23014LRa023L2_01.GSvar.gz -add 23014LRa023L2_01
	php ../compare_qcml_xml.php 23014LRa023L2_01_stats_fastq.qcML.out 23014LRa023L2_01_stats_fastq.qcML
	php ../compare_qcml_xml.php 23014LRa023L2_01_stats_map.qcML.out 23014LRa023L2_01_stats_map.qcML
	php ../compare_qcml_xml.php 23014LRa023L2_01_stats_vc.qcML.out 23014LRa023L2_01_stats_vc.qcML 1.01
	php ../compare_qcml_xml.php 23014LRa023L2_01_stats_other.qcML.out 23014LRa023L2_01_stats_other.qcML
	diff 23014LRa023L2_01_cnvs_clincnv.tsv.out 23014LRa023L2_01_cnvs_clincnv.tsv --ignore-matching-lines="#"
	php ../compare_tsv.php -in1 23014LRa023L2_01_var_structural_variants.bedpe.out -in2 23014LRa023L2_01_var_structural_variants.bedpe -skip_comments_matching command=,commandline=,fileDate=,PIPELINE=,SAMPLE= -skip_cols GENES,GENES_BREAKPOINTS,GENE_INFO,NGSD_HOM,NGSD_HET,NGSD_AF,NGSD_group,NGSD_SV_BREAKPOINT_DENSITY,NGSD_SV_BREAKPOINT_DENSITY_SYS,NGSD_PATHOGENIC_SVS
	diff 23014LRa023L2_01_prs.tsv.out 23014LRa023L2_01_prs.tsv
	diff 23014LRa023L2_01_phasing_track.bed.out 23014LRa023L2_01_phasing_track.bed
	diff 23014LRa023L2_01_bafs.igv.out 23014LRa023L2_01_bafs.igv
	diff 23014LRa023L2_01.paraphase.json.out paraphase/23014LRa023L2_01.paraphase.json

# 2nd test with BAM input (default test)
bam_analysis:
	php ../link_test_data.php 23014LRa023_01.mod.unmapped.bam . -copy
	php ../../src/Pipelines/analyze_longread.php -folder . -name 23014LRa023_01 -system system.ini -steps ma,vc,cn,sv,ph,re,me,pg,an -threads 4 -test --log analyze_longread_bam_all_$(TIMESTAMP).log

bam_diff:
	php ../compare_variants.php -in1 23014LRa023_01.GSvar.gz.out -in2 23014LRa023_01.GSvar.gz -add 23014LRa023_01
	php ../compare_qcml_xml.php 23014LRa023_01_stats_fastq.qcML.out 23014LRa023_01_stats_fastq.qcML
	php ../compare_qcml_xml.php 23014LRa023_01_stats_map.qcML.out 23014LRa023_01_stats_map.qcML
	php ../compare_qcml_xml.php 23014LRa023_01_stats_vc.qcML.out 23014LRa023_01_stats_vc.qcML 1.01
	php ../compare_qcml_xml.php 23014LRa023_01_stats_other.qcML.out 23014LRa023_01_stats_other.qcML
	diff 23014LRa023_01_cnvs_clincnv.tsv.out 23014LRa023_01_cnvs_clincnv.tsv --ignore-matching-lines="#"
	php ../compare_tsv.php -in1 23014LRa023_01_var_structural_variants.bedpe.out -in2 23014LRa023_01_var_structural_variants.bedpe -skip_comments_matching command=,commandline=,fileDate=,PIPELINE=,SAMPLE= -skip_cols GENES,GENES_BREAKPOINTS,GENE_INFO,NGSD_HOM,NGSD_HET,NGSD_AF,NGSD_group,NGSD_SV_BREAKPOINT_DENSITY,NGSD_SV_BREAKPOINT_DENSITY_SYS,NGSD_PATHOGENIC_SVS
	diff 23014LRa023_01_prs.tsv.out 23014LRa023_01_prs.tsv
	diff 23014LRa023_01_phasing_track.bed.out 23014LRa023_01_phasing_track.bed
	diff 23014LRa023_01_bafs.igv.out 23014LRa023_01_bafs.igv
	diff 23014LRa023_01_var_methylation.tsv.out 23014LRa023_01_var_methylation.tsv
	diff 23014LRa023_01.paraphase.json.out paraphase/23014LRa023_01.paraphase.json
	diff 23014LRa023_01_epigen.tsv.out 23014LRa023_01_epigen.tsv

# 3rd test with each step separate (not performed by default)
bam_analysis_single_steps:
	php ../link_test_data.php 23014LRa023_01.mod.unmapped.bam . -copy
	php ../../src/Pipelines/analyze_longread.php -folder . -name 23014LRa023_01 -system system.ini -steps ma -threads 4 --log analyze_longread_bam_ma_$(TIMESTAMP).log
	php ../../src/Pipelines/analyze_longread.php -folder . -name 23014LRa023_01 -system system.ini -steps vc -threads 4 --log analyze_longread_bam_vc_$(TIMESTAMP).log
	php ../../src/Pipelines/analyze_longread.php -folder . -name 23014LRa023_01 -system system.ini -steps cn -threads 4 --log analyze_longread_bam_cn_$(TIMESTAMP).log -test
	php ../../src/Pipelines/analyze_longread.php -folder . -name 23014LRa023_01 -system system.ini -steps sv -threads 4 --log analyze_longread_bam_sv_$(TIMESTAMP).log
	php ../../src/Pipelines/analyze_longread.php -folder . -name 23014LRa023_01 -system system.ini -steps ph -threads 4 --log analyze_longread_bam_ph_$(TIMESTAMP).log
	php ../../src/Pipelines/analyze_longread.php -folder . -name 23014LRa023_01 -system system.ini -steps re -threads 4 --log analyze_longread_bam_re_$(TIMESTAMP).log
	php ../../src/Pipelines/analyze_longread.php -folder . -name 23014LRa023_01 -system system.ini -steps pg -threads 4 --log analyze_longread_bam_pg_$(TIMESTAMP).log
	php ../../src/Pipelines/analyze_longread.php -folder . -name 23014LRa023_01 -system system.ini -steps an -threads 4 --log analyze_longread_bam_an_$(TIMESTAMP).log
	php ../../src/Pipelines/analyze_longread.php -folder . -name 23014LRa023_01 -system system.ini -steps me -threads 4 --log analyze_longread_bam_me_$(TIMESTAMP).log -test
	
	
	


