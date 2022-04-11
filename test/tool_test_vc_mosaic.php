<?php

require_once("framework.php");



$name = "vc_mosaic";
start_test($name);


$out_file1 = output_folder().$name."_out1.vcf";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in1.bam -vcf ".data_folder().$name."_in1.vcf -out $out_file1 -type WES -genes ".data_folder().$name."_in_genes.txt -debug --log ".output_folder().$name."_out1.log");
remove_lines_containing($out_file1, ["contig=", "fileDate=", "commandline="]);
check_file($out_file1, data_folder().$name."_out1.vcf");


end_test();

?>


