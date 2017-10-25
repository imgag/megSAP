<?php

require_once("framework.php");

$name = "vcf2gsvar";
start_test($name);

//standard
$out_file1 = output_folder().$name."_out1.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in1.vcf -build GRCh37 -out $out_file1 --log ".output_folder().$name."_out1.log");
check_file($out_file1, data_folder().$name."_out1.GSvar");

//empty input VCF
$out_file2 = output_folder().$name."_out2.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in_empty.vcf -build GRCh37 -out $out_file2 --log ".output_folder().$name."_out2.log");
check_file($out_file2, data_folder().$name."_out2.GSvar", true);

//input file without genotype information (e.g. for multi-sample analysis)
$out_file3 = output_folder().$name."_out3.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in2.vcf -build GRCh37 -out $out_file3 -multi --log ".output_folder().$name."_out3.log");
check_file($out_file3, data_folder().$name."_out3.GSvar", true);

//standard
$out_file4 = output_folder().$name."_out4.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in1.vcf -build hg19 -out $out_file4 --log ".output_folder().$name."_out4.log");
check_file($out_file4, data_folder().$name."_out4.GSvar");

end_test();

?>
