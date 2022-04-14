<?php

require_once("framework.php");



$name = "vc_mosaic";
start_test($name);

//reference result file "vc_mosaic_out1.vcf" still contains a false positive along with the true positives.
//Decided by whether the muatations were called in the single sample calling of NA12878 ("/mnt/storage2/projects/test/WES_Onko_test/Sample_NA12878x2_29/NA12878x2_29.GSvar")

//false positive:
//chr17	78444641	.	G	T	0	.	MQM=60;SAP=3;ABP=0	GT:DP:AO:GQ	0/0:215:4:0

//base test:
$out_file1 = output_folder().$name."_out1.vcf";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in1.bam -vcf ".data_folder().$name."_in1.vcf -out $out_file1 -type WES -genes ".data_folder().$name."_in_genes.txt --log ".output_folder().$name."_out1.log");
remove_lines_containing($out_file1, ["contig=", "fileDate=", "commandline="]);
check_file($out_file1, data_folder().$name."_out1.vcf");

//set min obs:
$out_file2 = output_folder().$name."_out2.vcf";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in1.bam -vcf ".data_folder().$name."_in1.vcf -out $out_file2 -type WES -genes ".data_folder().$name."_in_genes.txt -min_obs 3 -extend 0 --log ".output_folder().$name."_out2.log");
remove_lines_containing($out_file2, ["contig=", "fileDate=", "commandline="]);
check_file($out_file2, data_folder().$name."_out2.vcf");

//multithreaded: 2 threads (same result as base test)
$out_file3 = output_folder().$name."_out3.vcf";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in1.bam -vcf ".data_folder().$name."_in1.vcf -out $out_file3 -type WES -genes ".data_folder().$name."_in_genes.txt -threads 2 --log ".output_folder().$name."_out1.log");
remove_lines_containing($out_file3, ["contig=", "fileDate=", "commandline="]);
check_file($out_file3, data_folder().$name."_out1.vcf");

//multithreaded: 4 threads
$out_file4 = output_folder().$name."_out4.vcf";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in1.bam -vcf ".data_folder().$name."_in1.vcf -out $out_file4 -type WES -genes ".data_folder().$name."_in_genes.txt -min_obs 3 -extend 0 -threads 4 --log ".output_folder().$name."_out2.log");
remove_lines_containing($out_file4, ["contig=", "fileDate=", "commandline="]);
check_file($out_file4, data_folder().$name."_out2.vcf");

end_test();

?>


