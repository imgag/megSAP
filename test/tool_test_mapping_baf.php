<?php

require_once("framework.php");

$name = "mapping_baf";
start_test($name);

$t_bam = data_folder()."vc_strelka2_tu_in.bam";
$n_bam = data_folder()."vc_strelka2_no_in.bam";

$base_vcf = data_folder()."/".$name."_snps.vcf.gz";
$target = data_folder()."/".$name."_target.bed";

// test1: target, non-WGS, somatic
// expected output: variants in target regions, two samples/columns
$out1 = output_folder().$name."_out1.igv";
check_exec("php ".src_folder()."/NGS/mapping_baf.php -in $t_bam -n_in $n_bam -out $out1 -target $target -base_vcf $base_vcf --log ".output_folder().$name."_1.log");
check_file($out1, data_folder().$name."_out1.igv");

// test2: target, non-WGS
// expected output: variants in target regions, one sample/column
$out2 = output_folder().$name."_out2.igv";
check_exec("php ".src_folder()."/NGS/mapping_baf.php -in $t_bam -out $out2 -target $target -base_vcf $base_vcf --log ".output_folder().$name."_2.log");
check_file($out2, data_folder().$name."_out2.igv");

// test3: target, WGS
// expected output: 1st variant only
$out3 = output_folder().$name."_out3.igv";
check_exec("php ".src_folder()."/NGS/mapping_baf.php -in $t_bam -out $out3 -target $target -wgs -base_vcf $base_vcf --log ".output_folder().$name."_3.log");
check_file($out3, data_folder().$name."_out3.igv");

// test4: no target, WGS
// expected output: same results as test 3
$out4 = output_folder().$name."_out4.igv";
check_exec("php ".src_folder()."/NGS/mapping_baf.php -in $t_bam -out $out4 -wgs -base_vcf $base_vcf --log ".output_folder().$name."_4.log");
check_file($out4, data_folder().$name."_out3.igv");

end_test();

?>