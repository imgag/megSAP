<?php

require_once("framework.php");

$name = "filter_vcf";

start_test("filter_vcf");

$out_file1 = output_folder().$name."_out1.vcf";
check_exec("php ".src_folder()."/NGS/filter_vcf.php -in ".data_folder().$name."_in1.vcf -out $out_file1 -keep -type all -contamination 0.03 -roi ../data/enrichment/SeqCapEZv2_2013_02_19.bed");
check_file($out_file1, data_folder().$name."_out1.vcf");

$out_file2 = output_folder().$name."_out2.vcf";
check_exec("php ".src_folder()."/NGS/filter_vcf.php -in ".data_folder().$name."_in2.vcf -out $out_file2 -keep -type not-coding-splicing,synonymous,somatic-lq,off-target -contamination 0.03 -roi ../data/enrichment/ssHAEv5_2016_07_11.bed");
check_file($out_file2, data_folder().$name."_out2.vcf");

$out_file3 = output_folder().$name."_out3.vcf";
$return = check_exec("php ".src_folder()."/NGS/filter_vcf.php -in ".data_folder().$name."_in3.vcf -out $out_file3 -keep -type not-coding-splicing -roi ../data/enrichment/ssHAEv5_2016_07_11.bed",FALSE);
list($return,) = explode(" in ",$return[0]);
check($return,"ERROR: 'Found 'not-cod-spli' filter multiple times. Looks like this files was already annotated by filter_vcf.'");

$out_file4 = output_folder().$name."_out4.vcf";
check_exec("php ".src_folder()."/NGS/filter_vcf.php -in ".data_folder().$name."_in4.vcf -out $out_file4 -type not-coding-splicing,synonymous,somatic-lq,off-target -roi ../data/enrichment/SeqCapEZv2_2013_02_19.bed");
check_file($out_file4, data_folder().$name."_out4.vcf");

$out_file5 = output_folder().$name."_out5.vcf";
$return = check_exec("php ".src_folder()."/NGS/filter_vcf.php -in ".data_folder().$name."_in5.vcf -out $out_file5 -type not-coding-splicing,somatic-lq,off-target -ignore_filter -keep -roi ../data/enrichment/SeqCapEZv2_2013_02_19.bed",FALSE);
list($return,) = explode(" in ",$return[0]);
check($return,"ERROR: 'Found 'not-cod-spli, off-target' filter multiple times. Looks like this files was already annotated by filter_vcf.'");

end_test();
?>
