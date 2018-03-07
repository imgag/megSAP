
all:
	make clear link vc an cn diff
	
clear:
	rm -rf *.bam *.bai *.log *.GSvar *.vcf.gz *.vcf.gz.tbi *.qcML

link:
	php ../link_test_data.php GS140127_01.bam
	php ../link_test_data.php GS140127_01.bam.bai
	php ../link_test_data.php GS140549_01.bam
	php ../link_test_data.php GS140549_01.bam.bai
	php ../link_test_data.php GS140550_01.bam
	php ../link_test_data.php GS140550_01.bam.bai

vc:
	rm -rf multi_log1_vc.log
	php ../../src/Pipelines/multisample.php -bams GS140127_01.bam GS140549_01.bam GS140550_01.bam -status affected control control -system system.ini -out_folder . --log multi_log1_vc.log -steps vc

an:
	rm -rf multi_log2_an.log
	php ../../src/Pipelines/multisample.php -bams GS140127_01.bam GS140549_01.bam GS140550_01.bam -status affected control control -system system.ini -out_folder . --log multi_log2_an.log -steps an

cn:
	rm -rf multi_log3_cn.log
	php ../../src/Pipelines/multisample.php -bams GS140127_01.bam GS140549_01.bam GS140550_01.bam -status affected control control -system system.ini -out_folder . --log multi_log3_cn.log -steps cn

diff:
	php ../compare_variants.php multi.GSvar.out multi.GSvar
	diff multi_cnvs.tsv multi_cnvs.tsv.out
