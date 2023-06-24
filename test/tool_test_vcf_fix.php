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

//mosaic mode:
$out_file3 = output_folder().$name."_out3.vcf";
check_exec("cat ".data_folder().$name."_in1.vcf | php ".src_folder()."/NGS/{$name}.php --mosaic_mode  > $out_file3");
check_file($out_file3, data_folder().$name."_out3.vcf");

//longread mode
$out_file4 = output_folder().$name."_out4.vcf";
check_exec("cat ".data_folder().$name."_in3.vcf | php ".src_folder()."/NGS/{$name}.php --longread_mode  > $out_file4");
check_file($out_file4, data_folder().$name."_out4.vcf");

end_test();

?>
