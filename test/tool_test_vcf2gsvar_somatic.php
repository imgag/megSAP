<?php

require_once("framework.php");

$name = "vcf2gsvar_somatic";
start_test($name);

//strelka (tumor-normal pair)
$out_file1 = output_folder().$name."_out1.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in1.vcf -out $out_file1 -t_col GS110168 -n_col GS110169 --log ".output_folder().$name."_out1.log");
check_file($out_file1, data_folder().$name."_out1.GSvar");

//freebayes (tumor-only)
$out_file2 = output_folder().$name."_out2.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in2.vcf -out $out_file2 -t_col GS150344_01 --log ".output_folder().$name."_out2.log");
check_file($out_file2, data_folder().$name."_out2.GSvar");

//freebayes (tumor-only) with refseq annotation
$out_file3 = output_folder().$name."_out3.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in3.vcf -out $out_file3 -t_col GS150344_01 --log ".output_folder().$name."_out3.log");
check_file($out_file3, data_folder().$name."_out3.GSvar");

end_test();

?>
