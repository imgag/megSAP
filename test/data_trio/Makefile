
all:
	make clear link vc an cn diff
	
clear:
	rm -rf *.bam *.bai *.log *.GSvar *.tsv *.vcf.gz *.vcf.gz.tbi *.qcML

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
	#all variants
	php ../compare_variants.php -in1 trio.GSvar.out -in2 trio.GSvar
	#variants after filtering (only when NGSD is enabled)
	php ../filter_variants.php -in trio.GSvar -out trio_denovo.GSvar -filters filter_denovo.txt -requires_db NGSD
	php ../compare_variants.php -in1 trio_denovo.GSvar.out -in2 trio_denovo.GSvar -requires_db NGSD
	php ../filter_variants.php -in trio.GSvar -out trio_recessive.GSvar -filters filter_recessive.txt -requires_db NGSD
	php ../compare_variants.php -in1 trio_recessive.GSvar.out -in2 trio_recessive.GSvar -requires_db NGSD
	php ../filter_variants.php -in trio.GSvar -out trio_comphet.GSvar -filters filter_comphet.txt -requires_db NGSD
	php ../compare_variants.php -in1 trio_comphet.GSvar.out -in2 trio_comphet.GSvar -requires_db NGSD
	php ../filter_variants.php -in trio.GSvar -out trio_loh.GSvar -filters filter_loh.txt -requires_db NGSD
	php ../compare_variants.php -in1 trio_loh.GSvar.out -in2 trio_loh.GSvar -requires_db NGSD
	#UPDs
	diff -b trio_upd.tsv.out trio_upd.tsv

### debug targets

quick_test:
	rm -rf  trio_quick.log
	php ../../src/Pipelines/trio.php -c GS140127_01.bam -m GS140549_01.bam -f GS140550_01.bam -out_folder . --log trio_quick.log -no_check -system system_quick.ini
