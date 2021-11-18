<?php

require_once("framework.php");

$name = "vcf2gsvar";
start_test($name);

//genotype_mode=single
$out_file1 = output_folder().$name."_out1.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in1.vcf -out $out_file1 --log ".output_folder().$name."_out1.log");
check_file($out_file1, data_folder().$name."_out1.GSvar");

//genotype_mode=skip, with upstream/downstream 
$out_file2 = output_folder().$name."_out2.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in1.vcf -updown -out $out_file2 -genotype_mode skip --log ".output_folder().$name."_out2.log");
check_file($out_file2, data_folder().$name."_out2.GSvar");

//genotype_mode=single, NGSD group counts
$out_file_db = output_folder().$name."_out_db2.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in_db2.vcf -out $out_file_db --log ".output_folder().$name."_out_db2.log");
check_file($out_file_db, data_folder().$name."_out_db2.GSvar");

//genotype_mode=single, empty input VCF
$out_file_empty = output_folder().$name."_out_empty.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in_empty.vcf -out $out_file_empty --log ".output_folder().$name."_out_empty.log");
remove_lines_containing($out_file_empty, "#CREATION_DATE=");
check_file($out_file_empty, data_folder().$name."_out_empty.GSvar", true);

//genotype_mode=single, WGS mode
$out_file1 = output_folder().$name."_out5.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in1.vcf -updown -wgs -out $out_file1 --log ".output_folder().$name."_out1.log");
check_file($out_file1, data_folder().$name."_out5.GSvar");

//genotype_mode=multi
$out_file_multi = output_folder().$name."_out_multi.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in_multi.vcf -out $out_file_multi -genotype_mode multi --log ".output_folder().$name."_out_multi.log");
remove_lines_containing($out_file_multi, "#CREATION_DATE=");
check_file($out_file_multi, data_folder().$name."_out_multi.GSvar", true);

//genotype_mode=single, DRAGEN
$out_file_dragen = output_folder().$name."_out_dragen.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in_dragen.vcf -out $out_file_dragen --log ".output_folder().$name."_out_dragen.log");
check_file($out_file_dragen, data_folder().$name."_out_dragen.GSvar");

end_test();

?>
