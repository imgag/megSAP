<?php

require_once("framework.php");

$name = "vcf2gsvar";
start_test($name);

//genotype_mode=single
$out_file1 = output_folder().$name."_out1.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in1.vcf -blacklist -out $out_file1 --log ".output_folder().$name."_out1.log");
check_file($out_file1, data_folder().$name."_out1.GSvar");

//genotype_mode=single, empty input VCF
$out_file2 = output_folder().$name."_out2.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in_empty.vcf -blacklist -out $out_file2 --log ".output_folder().$name."_out2.log");
check_file($out_file2, data_folder().$name."_out2.GSvar", true);

//genotype_mode=multi
$out_file3 = output_folder().$name."_out3.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in2.vcf -blacklist -out $out_file3 -genotype_mode multi --log ".output_folder().$name."_out3.log");
check_file($out_file3, data_folder().$name."_out3.GSvar", true);

//genotype_mode=skip, with upstream/downstream 
$out_file4 = output_folder().$name."_out4.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in1.vcf -blacklist  -updown -out $out_file4 -genotype_mode skip --log ".output_folder().$name."_out4.log");
check_file($out_file4, data_folder().$name."_out4.GSvar");

//genotype_mode=single, no blacklist
$out_file5 = output_folder().$name."_out5.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in1.vcf -out $out_file5 --log ".output_folder().$name."_out5.log");
check_file($out_file5, data_folder().$name."_out5.GSvar");

//genotype_mode=single, DRAGEN
$out_file6 = output_folder().$name."_out6.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in_dragen.vcf -out $out_file6 --log ".output_folder().$name."_out6.log");
check_file($out_file6, data_folder().$name."_out6.GSvar");

//genotype_mode=single, RefSeq annotation
$out_file7 = output_folder().$name."_out7.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in3.vcf -blacklist -out $out_file7 --log ".output_folder().$name."_out7.log");
check_file($out_file7, data_folder().$name."_out7.GSvar");

//genotype_mode=single, RefSeq annotation
$out_file8 = output_folder().$name."_out8.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in4.vcf -blacklist -out $out_file8 --log ".output_folder().$name."_out8.log");
check_file($out_file8, data_folder().$name."_out8.GSvar");


end_test();

?>
