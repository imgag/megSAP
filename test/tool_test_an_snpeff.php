<?php

require_once("framework.php");

$name = "an_snpeff";
start_test($name);

//standard
$out_file1 = output_folder().$name."_out1.vcf";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in1.vcf -out $out_file1 --log ".output_folder().$name."_out1.log");
remove_lines_containing($out_file1, array("SnpSiftCmd"));
check_file($out_file1, data_folder().$name."_out1.vcf", true);

//empty input VCF
$out_file2 = output_folder().$name."_out2.vcf";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in_empty.vcf -out $out_file2 --log ".output_folder().$name."_out2.log");
remove_lines_containing($out_file2, array("SnpSiftCmd"));
check_file($out_file2, data_folder().$name."_out2.vcf", true);

//standard - no upstream/downstream
$out_file3 = output_folder().$name."_out3.vcf";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in1.vcf -out $out_file3 --log ".output_folder().$name."_out3.log -no_updown");
remove_lines_containing($out_file3, array("SnpSiftCmd"));
check_file($out_file3, data_folder().$name."_out3.vcf", true);

end_test();

?>
