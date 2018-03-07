
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
	php ../compare_variants.php trio.GSvar.out trio.GSvar
	cut -f 1-5,9 trio.GSvar | egrep -v "^#|off-target" | egrep "trio_" > trio_special.tsv
	diff -b trio_special.tsv.out trio_special.tsv

### debug targets

quick_test:
	rm -rf  trio_quick.log
	php ../../src/Pipelines/trio.php -c GS140127_01.bam -m GS140549_01.bam -f GS140550_01.bam -out_folder . --log trio_quick.log -no_check -system system_quick.ini
