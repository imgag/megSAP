<?php

require_once("framework.php");

$name = "vcf2gsvar_somatic";
start_test($name);

//init
check_exec(get_path("ngs-bits")."NGSDInit -test -add ".data_folder()."/{$name}.sql");


//strelka (tumor-normal pair)
$out_file1 = output_folder().$name."_out1.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in1.vcf -out $out_file1 -t_col GS110168_03 -n_col GS110169_03 -db NGSD_TEST --log ".output_folder().$name."_out1.log");
check_file($out_file1, data_folder().$name."_out1.GSvar");

//freebayes (tumor-only)
$out_file2 = output_folder().$name."_out2.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in2.vcf -out $out_file2 -t_col GS150344_01 -db NGSD_TEST --log ".output_folder().$name."_out2.log");
check_file($out_file2, data_folder().$name."_out2.GSvar");

//freebayes (tumor-only) with refseq annotation
$out_file3 = output_folder().$name."_out3.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in3.vcf -out $out_file3 -t_col GS150344_01 -db NGSD_TEST --log ".output_folder().$name."_out3.log");
check_file($out_file3, data_folder().$name."_out3.GSvar");

//varscan2 (tumor_only)
$out_file4 = output_folder().$name."_out4.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in4.vcf -out $out_file4 -t_col GS150344_01 -db NGSD_TEST --log ".output_folder().$name."_out4.log");
check_file($out_file4, data_folder().$name."_out4.GSvar");

//varscan2 (tumor only) and COSMIC annotation
$out_file5 = output_folder().$name."_out5.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in5.vcf -out $out_file5 -t_col GS150344_01 -db NGSD_TEST --log ".output_folder().$name."_out5.log");
check_file($out_file5, data_folder().$name."_out5.GSvar");

end_test();

?>
