<?php

require_once("framework.php");

$name = "mapping_baf";
start_test($name);

$t_bam = data_folder()."vc_strelka2_tu_in.bam";
$n_bam = data_folder()."vc_strelka2_no_in.bam";

$base_vcf = data_folder()."/".$name."_snps.vcf.gz";
$target = data_folder()."/".$name."_target.bed";

// test1: somatic, target 
// expected output: variants in target regions, two samples/columns
$out1 = output_folder().$name."_out1.igv";
check_exec("php ".src_folder()."/NGS/mapping_baf.php -in $t_bam -n_in $n_bam -out $out1 -target $target -base_vcf $base_vcf --log ".output_folder().$name."_1.log");
check_file($out1, data_folder().$name."_out1.igv");

// test2: germline, target 
// expected output: variants in target regions, one sample/column
$out2 = output_folder().$name."_out2.igv";
check_exec("php ".src_folder()."/NGS/mapping_baf.php -in $t_bam -out $out2 -target $target -base_vcf $base_vcf --log ".output_folder().$name."_2.log");
check_file($out2, data_folder().$name."_out2.igv");

// test3: germline, no target => WGS
// expected output: 1st variant only
$out3 = output_folder().$name."_out3.igv";
check_exec("php ".src_folder()."/NGS/mapping_baf.php -in $t_bam -out $out3 -base_vcf $base_vcf --log ".output_folder().$name."_3.log");
check_file($out3, data_folder().$name."_out3.igv");


end_test();

?>