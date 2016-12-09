<?php

require_once("framework.php");

$name = "vcf_fix";
start_test($name);

//single-sample input
$out_file1 = output_folder().$name."_out1.vcf";
check_exec("cat ".data_folder().$name."_in1.vcf | php ".src_folder()."/NGS/{$name}.php > $out_file1");
check_file($out_file1, data_folder().$name."_out1.vcf");

//multi-sample input
$out_file2 = output_folder().$name."_out2.vcf";
check_exec("cat ".data_folder().$name."_in2.vcf | php ".src_folder()."/NGS/{$name}.php > $out_file2");
check_file($out_file2, data_folder().$name."_out2.vcf");

end_test();

?>
