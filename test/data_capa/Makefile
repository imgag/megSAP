all: clean link complete_pipeline report_only diff

clean:
	rm -f *.bam *.bai *.log *.txt *.GSvar *.vcf *.fastq.gz *.qcML *_qcML.txt *_igv.xml *.log *_report.rtf *.bed
	rm -f Somatic_Tumor1-Normal1/Tumor1-Normal1.log Somatic_Tumor1-Normal1/Tumor1-Normal1_igv.xml Somatic_Tumor1-Normal1/Tumor1-Normal1_report.txt
	rm -f Somatic_Tumor2-Normal2/Tumor2-Normal2.log Somatic_Tumor2-Normal2/Tumor2-Normal2_igv.xml Somatic_Tumor2-Normal2/Tumor2-Normal2_report.txt
	rm -rf Sample_GS150496_01
	rm -rf Sample_GS150497_01
	rm -rf Somatic_GS150497_01-GS150496_01

link:
	mkdir -p Sample_GS150496_01
	php ../link_test_data.php GS150496_01_S1_L001extracted_R1_001.fastq.gz Sample_GS150496_01/GS150496_01_S1_L001extracted_R1_001.fastq.gz
	php ../link_test_data.php GS150496_01_S1_L001extracted_R2_001.fastq.gz Sample_GS150496_01/GS150496_01_S1_L001extracted_R2_001.fastq.gz
	mkdir -p Sample_GS150497_01
	php ../link_test_data.php GS150497_01_S2_L001extracted_R1_001.fastq.gz Sample_GS150497_01/GS150497_01_S2_L001extracted_R1_001.fastq.gz
	php ../link_test_data.php GS150497_01_S2_L001extracted_R2_001.fastq.gz Sample_GS150497_01/GS150497_01_S2_L001extracted_R2_001.fastq.gz

complete_pipeline:
	mkdir -p Somatic_GS150497_01-GS150496_01
	php ../../src/Pipelines/somatic_capa.php -t_id GS150497_01 -t_sys system.ini -n_id GS150496_01 -n_sys system.ini -p_folder . -o_folder Somatic_GS150497_01-GS150496_01 -steps ma,vc,an,ci,re -nsc --log GS150497_01-GS150496_01.log
	cp Somatic_GS150497_01-GS150496_01/GS150497_01-GS150496_01.GSvar GS150497_01-GS150496_01.GSvar
	cp Somatic_GS150497_01-GS150496_01/GS150497_01-GS150496_01_report.rtf GS150497_01-GS150496_01_report.rtf
	cp Somatic_GS150497_01-GS150496_01/GS150497_01-GS150496_01_stat_lowcov.bed GS150497_01-GS150496_01_stat_lowcov.bed

report_only:
	php ../../src/Pipelines/somatic_capa.php -t_id Tumor1 -t_sys system.ini -n_id Normal1 -n_sys system.ini -p_folder . -o_folder Somatic_Tumor1-Normal1 -steps re --log Somatic_Tumor1-Normal1/Tumor1-Normal1.log
	php ../../src/Pipelines/somatic_capa.php -t_id Tumor2 -t_sys system.ini -n_id Normal2 -n_sys system.ini -p_folder . -o_folder Somatic_Tumor2-Normal2 -steps re --log Somatic_Tumor2-Normal2/Tumor2-Normal2.log
	cp Somatic_Tumor1-Normal1/Tumor1-Normal1_report.rtf Tumor1-Normal1_report.rtf
	cp Somatic_Tumor1-Normal1/Tumor1-Normal1_igv.xml Tumor1-Normal1_igv.xml
	cp Somatic_Tumor2-Normal2/Tumor2-Normal2_report.rtf Tumor2-Normal2_report.rtf
	cp Somatic_Tumor2-Normal2/Tumor2-Normal2_igv.xml Tumor2-Normal2_igv.xml

diff:
	php ../compare_variants.php GS150497_01-GS150496_01.GSvar.out GS150497_01-GS150496_01.GSvar
	diff -b -I ' Wiss. Bericht GS150497_01-GS150496_01 vom ' -I 'Datum:' -I 'Revision der Analysepipeline' -I 'Tumoranteil histol./molekular:' GS150497_01-GS150496_01_report.rtf.out GS150497_01-GS150496_01_report.rtf
	ls */GS*.qcML | sort | xargs grep "<qualityParameter" > somatic_qcML.txt
	php ../compare_qcml.php somatic_qcML.txt.out somatic_qcML.txt
	diff -b -I ' Wiss. Bericht Tumor1-Normal1 vom ' -I 'Datum:' -I 'Revision der Analysepipeline' -I 'Tumoranteil histol./molekular:' Tumor1-Normal1_report.rtf.out Tumor1-Normal1_report.rtf
	diff -b -I ' Wiss. Bericht Tumor2-Normal2 vom ' -I 'Datum:' -I 'Revision der Analysepipeline' -I 'Tumoranteil histol./molekular:' Tumor2-Normal2_report.rtf.out Tumor2-Normal2_report.rtf
	diff -b Tumor1-Normal1_igv.xml.out Tumor1-Normal1_igv.xml
	diff -b Tumor2-Normal2_igv.xml.out Tumor2-Normal2_igv.xml
	diff -b -I '#ROI:' GS150497_01-GS150496_01_stat_lowcov.bed.out GS150497_01-GS150496_01_stat_lowcov.bed
