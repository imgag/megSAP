<?php

require_once("framework.php");

$name = "vc_viral_load";
start_test($name);

$out_vi_bam = output_folder().$name."_out1.bam";
$out_vi_tsv = output_folder().$name."_out1.tsv";
check_exec("php ".src_folder()."/Tools/{$name}.php -in ".data_folder().$name."_in1.bam -avg_target_cov 180.11 -viral_bam $out_vi_bam -viral_cov $out_vi_tsv --log ".output_folder().$name."_out1.log");
check_file($out_vi_tsv, data_folder().$name."_out1.tsv", true);



end_test();

?>
