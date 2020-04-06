<?php

require_once("framework.php");

$name = "vcf2gsvar";
start_test($name);

//legacy test
$out_file1 = output_folder().$name."_out1_old.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in1_old.vcf -blacklist -out $out_file1 --log ".output_folder().$name."_out1_old.log");
check_file($out_file1, data_folder().$name."_out1_old.GSvar");

//genotype_mode=single
$out_file1 = output_folder().$name."_out1.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in1.vcf -blacklist -out $out_file1 --log ".output_folder().$name."_out1.log");
check_file($out_file1, data_folder().$name."_out1.GSvar");

//genotype_mode=skip, with upstream/downstream 
$out_file2 = output_folder().$name."_out2.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in1.vcf -blacklist  -updown -out $out_file2 -genotype_mode skip --log ".output_folder().$name."_out2.log");
check_file($out_file2, data_folder().$name."_out2.GSvar");

//genotype_mode=single, no blacklist
$out_file3 = output_folder().$name."_out3.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in1.vcf -out $out_file3 --log ".output_folder().$name."_out3.log");
check_file($out_file3, data_folder().$name."_out3.GSvar");

//genotype_mode=single, RefSeq annotation
$out_file4 = output_folder().$name."_out4.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in2.vcf -blacklist -out $out_file4 --log ".output_folder().$name."_out4.log");
check_file($out_file4, data_folder().$name."_out4.GSvar");

//genotype_mode=single, RefSeq annotation, NGSD group counts
$out_file_db = output_folder().$name."_out_db2.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in_db2.vcf -blacklist -out $out_file_db --log ".output_folder().$name."_out_db2.log");
check_file($out_file_db, data_folder().$name."_out_db2.GSvar");

//genotype_mode=single, empty input VCF
$out_file_empty = output_folder().$name."_out_empty.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in_empty.vcf -blacklist -out $out_file_empty --log ".output_folder().$name."_out_empty.log");
check_file($out_file_empty, data_folder().$name."_out_empty.GSvar", true);

//genotype_mode=multi
$out_file_multi = output_folder().$name."_out_multi.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in_multi.vcf -blacklist -out $out_file_multi -genotype_mode multi --log ".output_folder().$name."_out_multi.log");
check_file($out_file_multi, data_folder().$name."_out_multi.GSvar", true);

//genotype_mode=single, DRAGEN
$out_file_dragen = output_folder().$name."_out_dragen.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in_dragen.vcf -out $out_file_dragen --log ".output_folder().$name."_out_dragen.log");
check_file($out_file_dragen, data_folder().$name."_out_dragen.GSvar");

end_test();

?>
