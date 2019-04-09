
all:
	make clear link vc an cn diff
	
clear:
	rm -rf *.bam *.bai *.log *.GSvar *.vcf.gz *.vcf.gz.tbi *.qcML trio_cnvs_clincnv.tsv trio_upd.tsv

link:
	php ../link_test_data.php GS140127_01.bam
	php ../link_test_data.php GS140127_01.bam.bai
	php ../link_test_data.php GS140549_01.bam
	php ../link_test_data.php GS140549_01.bam.bai
	php ../link_test_data.php GS140550_01.bam
	php ../link_test_data.php GS140550_01.bam.bai

vc:
	rm -rf trio_vc.log
	php ../../src/Pipelines/trio.php -c GS140127_01.bam -m GS140549_01.bam -f GS140550_01.bam -system system.ini -out_folder . --log trio_vc.log -steps vc

an:
	rm -rf trio_an.log
	php ../../src/Pipelines/trio.php -c GS140127_01.bam -m GS140549_01.bam -f GS140550_01.bam -system system.ini -out_folder . --log trio_an.log -steps an

cn:
	rm -rf trio_cn.log
	php ../../src/Pipelines/trio.php -c GS140127_01.bam -m GS140549_01.bam -f GS140550_01.bam -system system.ini -out_folder . --log trio_cn.log -steps cn

diff:
	#variants
	php ../compare_variants.php -in1 trio.GSvar.out -in2 trio.GSvar
	#CNVs
	diff trio_cnvs_clincnv.tsv trio_cnvs_clincnv.tsv.out
	#UPDs
	diff -b trio_upd.tsv.out trio_upd.tsv
	#variants after filtering (only when NGSD is enabled)
	php ../filter_variants.php -in trio.GSvar -out trio_denovo.GSvar -filters filter_denovo.txt -requires_db NGSD
	php ../compare_variants.php -in1 trio_denovo.GSvar.out -in2 trio_denovo.GSvar -requires_db NGSD
	
######################### AUXILARY TARGETS #########################

update_cnv_data:
	php ../../src/Pipelines/analyze.php -steps cn -name GS140127_01 -folder /mnt/projects/diagnostic/KM/Sample_GS140127_01/
	php ../../src/Pipelines/analyze.php -steps cn -name GS140549_01 -folder /mnt/projects/diagnostic/KM/Sample_GS140549_01/
	php ../../src/Pipelines/analyze.php -steps cn -name GS140550_01 -folder /mnt/projects/diagnostic/KM/Sample_GS140550_01/
	cp /mnt/projects/diagnostic/KM/Sample_GS140127_01/GS140127_01_cnvs_clincnv.tsv .
	cp /mnt/projects/diagnostic/KM/Sample_GS140549_01/GS140549_01_cnvs_clincnv.tsv .
	cp /mnt/projects/diagnostic/KM/Sample_GS140550_01/GS140550_01_cnvs_clincnv.tsv .
	cp /mnt/projects/diagnostic/KM/Sample_GS140127_01/GS140127_01_cnvs_clincnv.tsv ../data_multi/
	cp /mnt/projects/diagnostic/KM/Sample_GS140549_01/GS140549_01_cnvs_clincnv.tsv ../data_multi/
	cp /mnt/projects/diagnostic/KM/Sample_GS140550_01/GS140550_01_cnvs_clincnv.tsv ../data_multi/
