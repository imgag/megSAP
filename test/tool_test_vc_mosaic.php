<?php

require_once("framework.php");



$name = "vc_mosaic";
start_test($name);

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



end_test();

?>


