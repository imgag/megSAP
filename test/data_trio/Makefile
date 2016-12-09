
all:
	make clear link analyze diff
	
clear:
	rm -rf *.bam *.bai *.log *.GSvar *.tsv *.vcf.gz *.vcf.gz.tbi *.qcML

link:
	php ../link_test_data.php GS140127_01.bam
	php ../link_test_data.php GS140127_01.bam.bai
	php ../link_test_data.php GS140549_01.bam
	php ../link_test_data.php GS140549_01.bam.bai
	php ../link_test_data.php GS140550_01.bam
	php ../link_test_data.php GS140550_01.bam.bai

analyze:
	php ../../src/Pipelines/trio.php -c GS140127_01.bam -m GS140549_01.bam -f GS140550_01.bam -system system.ini -out_folder . --log trio.log

diff:
	php ../compare_variants.php GS140127_01.GSvar.out GS140127_01.GSvar
	cut -f 1,2,3,4,5,7 GS140127_01.GSvar | egrep -v "^#" | egrep "trio_" > GS140127_trio.tsv
	diff -b GS140127_trio.tsv.out GS140127_trio.tsv

### debug targets

analyze_start_at_vc:
	rm -rf trio.log
	php ../../src/Pipelines/trio.php -c GS140127_01.bam -m GS140549_01.bam -f GS140550_01.bam -system system.ini -out_folder . --log trio.log -start vc

analyze_start_at_an:
	rm -rf trio.log
	php ../../src/Pipelines/trio.php -c GS140127_01.bam -m GS140549_01.bam -f GS140550_01.bam -system system.ini -out_folder . --log trio.log -start an

quick_test:
	rm -rf  trio_quick.log
	php ../../src/Pipelines/trio.php -c GS140127_01.bam -m GS140549_01.bam -f GS140550_01.bam -out_folder . --log trio_quick.log -start realign -system system_quick.ini

quick_test_norealign:
	rm -rf  trio_quick_norealign.log
	php ../../src/Pipelines/trio.php -c GS140127_01.bam -m GS140549_01.bam -f GS140550_01.bam -out_folder . --log trio_quick_norealign.log -start vc -system system_quick.ini
