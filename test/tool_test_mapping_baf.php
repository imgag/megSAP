<?php

require_once("framework.php");

$name = "mapping_baf";
start_test($name);

$t_bam = data_folder()."vc_strelka_tu_in.bam";
$n_bam = data_folder()."vc_strelka_no_in.bam";

$out_file = output_folder().$name."_out1.igv";
check_exec("php ".src_folder()."/NGS/mapping_baf.php -in ".data_folder()."/vc_strelka2_tu_in.bam -out $out_file -test ".data_folder()."/mapping_baf_snps.vcf --log ".output_folder().$name."_1.log");
check_file($out_file, data_folder().$name."_out1.igv");

$out_file = output_folder().$name."_out2.igv";
check_exec("php ".src_folder()."/NGS/mapping_baf.php -in ".data_folder()."/vc_strelka2_tu_in.bam -n_in ".data_folder()."/vc_strelka2_no_in.bam -out $out_file  -test ".data_folder()."/mapping_baf_snps.vcf --log ".output_folder().$name."_2.log");
check_file($out_file, data_folder().$name."_out2.igv");

end_test();
?>
