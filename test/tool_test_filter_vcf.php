<?php

require_once("framework.php");

$name = "filter_vcf";

start_test("filter_vcf");

$out_file1 = output_folder().$name."_out1.vcf";
check_exec("php ".src_folder()."/NGS/filter_vcf.php -in ".data_folder().$name."_in1.vcf -out $out_file1 -keep -type all -contamination 0.03 -roi ../data/enrichment/SeqCapEZv2_2013_02_19.bed");
check_file($out_file1, data_folder().$name."_out1.vcf");

$out_file2 = output_folder().$name."_out2.vcf";
check_exec("php ".src_folder()."/NGS/filter_vcf.php -in ".data_folder().$name."_in2.vcf -out $out_file2 -keep -type iVac -contamination 0.03 -roi ../data/enrichment/ssHAEv5_2016_07_11.bed");
check_file($out_file2, data_folder().$name."_out2.vcf");

$out_file3 = output_folder().$name."_out3.vcf";
check_exec("php ".src_folder()."/NGS/filter_vcf.php -in ".data_folder().$name."_in3.vcf -out $out_file3 -keep -type coding -roi ../data/enrichment/ssHAEv5_2016_07_11.bed");
check_file($out_file3, data_folder().$name."_out3.vcf");

$out_file4 = output_folder().$name."_out4.vcf";
check_exec("php ".src_folder()."/NGS/filter_vcf.php -in ".data_folder().$name."_in4.vcf -out $out_file4 -type iVac -roi ../data/enrichment/SeqCapEZv2_2013_02_19.bed");
check_file($out_file4, data_folder().$name."_out4.vcf");

end_test();
?>
