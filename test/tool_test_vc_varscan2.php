<?php

require_once("framework.php");

$name = "vc_varscan2";
start_test($name);

$bam = data_folder().$name."_tu.bam";

$out_file = output_folder().$name."_1_out.vcf.gz";
check_exec("php ".src_folder()."/NGS/vc_varscan2.php -bam $bam -out $out_file -target ". data_folder().$name .".bed -name varscan2_tumor --log ".output_folder().$name."_1.log");
check_file($out_file, data_folder().$name."_1_out.vcf.gz");


end_test();
?>
