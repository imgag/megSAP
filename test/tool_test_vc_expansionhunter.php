<?php

require_once("framework.php");

$name = "vc_expansionhunter";
start_test($name);

//test
$out_file_vcf = output_folder().$name."_out1.vcf";
$out_file_json = output_folder().$name."_out1.json";
check_exec("php ".src_folder()."/NGS/".$name.".php -in ".data_folder().$name."_in1.bam -out $out_file_vcf");
remove_lines_containing($out_file_vcf, array("##filedate="));
check_file($out_file_vcf, data_folder().$name."_out1.vcf");
check_file($out_file_json, data_folder().$name."_out1.json");

end_test();

?>