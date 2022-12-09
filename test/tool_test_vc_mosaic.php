<?php

require_once("framework.php");



$name = "vc_mosaic";
start_test($name);

//reference result file "vc_mosaic_out1.vcf" still contains a false positive along with the true positives.
//Decided by whether the muatations were called in the single sample calling of NA12878 ("/mnt/storage2/projects/test/WES_Onko_test/Sample_NA12878x2_29/NA12878x2_29.GSvar")
//false positive:
//chr17	78444641	.	G	T	0	.	MQM=60;SAP=3;ABP=0	GT:DP:AO:GQ	0/0:215:4:0

//base test
$out_file1 = output_folder().$name."_out1.vcf.gz";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in1.bam -out $out_file1 -min_obs 2 -target ".data_folder().$name."_in_genes_ext_merged.bed --log ".output_folder().$name."_out1.log");
remove_lines_containing($out_file1, ["contig=", "fileDate=", "commandline="]);
check_file($out_file1, data_folder().$name."_out1.vcf.gz");

//set min obs and without zipping
$out_file2 = output_folder().$name."_out2.vcf.gz";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in1.bam -out $out_file2 -min_obs 3 -no_zip -target ".data_folder().$name."_in_genes_ext_merged.bed --log ".output_folder().$name."_out2.log");
remove_lines_containing($out_file2, ["contig=", "fileDate=", "commandline="]);
check_file($out_file2, data_folder().$name."_out2.vcf");

//multithreaded: 2 threads
$out_file3 = output_folder().$name."_out3.vcf.gz";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in1.bam -out $out_file3 -min_obs 2 -target ".data_folder().$name."_in_genes_ext_merged.bed -threads 2 --log ".output_folder().$name."_out3.log");
remove_lines_containing($out_file3, ["contig=", "fileDate=", "commandline="]);
check_file($out_file3, data_folder().$name."_out1.vcf.gz");

//multithreaded: 4 threads
$out_file4 = output_folder().$name."_out4.vcf.gz";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in1.bam -out $out_file4 -min_obs 2 -target ".data_folder().$name."_in_genes_ext_merged.bed  -threads 4 --log ".output_folder().$name."_out4.log");
remove_lines_containing($out_file4, ["contig=", "fileDate=", "commandline="]);
check_file($out_file4, data_folder().$name."_out1.vcf.gz");

end_test();

?>


