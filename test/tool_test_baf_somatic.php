<?php

require_once("framework.php");

$name = "baf_somatic";
start_test($name);

$t_bam = data_folder()."vc_strelka2_tu_in.bam";
$n_bam = data_folder()."vc_strelka2_no_in.bam";
$vcf = data_folder()."/".$name."_sites.vcf.gz";

//test 1
$out1 = output_folder().$name."_out1.igv";
$log1 = output_folder().$name."_out1.log";
check_exec("php ".src_folder()."/NGS/baf_somatic.php -bam_t {$t_bam} -bam_n {$n_bam} -vcf {$vcf} -out {$out1} -depth --log {$log1}");
check_file($out1, data_folder().$name."_out1.igv");

//test 2 - downsampling
$out2 = output_folder().$name."_out2.igv";
$log2 = output_folder().$name."_out2.log";
check_exec("php ".src_folder()."/NGS/baf_somatic.php -bam_t {$t_bam} -bam_n {$n_bam} -vcf {$vcf} -out {$out2} -depth -downsample 2 --log {$log2}");
check_file($out2, data_folder().$name."_out2.igv");

$out3 = output_folder().$name."_out3.igv";
$log3 = output_folder().$name."_out3.log";
check_exec("php ".src_folder()."/NGS/baf_somatic.php -bam_t {$t_bam} -bam_n {$n_bam} -vcf {$vcf} -out {$out3} -depth -downsample 3 --log {$log3}");
check_file($out3, data_folder().$name."_out3.igv");

end_test();
?>