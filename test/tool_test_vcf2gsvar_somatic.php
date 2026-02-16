<?php

require_once("framework.php");

$name = "vcf2gsvar_somatic";
start_test($name);

//strelka (tumor-normal pair), vcf file has COSMIC CMC and CANCERHOTSPOTS annotation
$out_file1 = output_folder().$name."_out1.GSvar";
check_exec("php ".src_folder()."/Tools/{$name}.php -in ".data_folder().$name."_in1.vcf -out $out_file1 -t_col GS110168_03 -n_col GS110169_03 -db NGSD_TEST --log ".output_folder().$name."_out1.log");
check_file($out_file1, data_folder().$name."_out1.GSvar");

//freebayes (tumor-only)
$out_file2 = output_folder().$name."_out2.GSvar";
check_exec("php ".src_folder()."/Tools/{$name}.php -in ".data_folder().$name."_in2.vcf -out $out_file2 -t_col GS150344_01 -db NGSD_TEST --log ".output_folder().$name."_out2.log");
check_file($out_file2, data_folder().$name."_out2.GSvar");

//freebayes (tumor-only) with refseq annotation
$out_file3 = output_folder().$name."_out3.GSvar";
check_exec("php ".src_folder()."/Tools/{$name}.php -in ".data_folder().$name."_in3.vcf -out $out_file3 -t_col GS150344_01 -db NGSD_TEST --log ".output_folder().$name."_out3.log");
check_file($out_file3, data_folder().$name."_out3.GSvar");

//varscan2 (tumor_only)
$out_file4 = output_folder().$name."_out4.GSvar";
check_exec("php ".src_folder()."/Tools/{$name}.php -in ".data_folder().$name."_in4.vcf -out $out_file4 -t_col GS150344_01 -db NGSD_TEST --log ".output_folder().$name."_out4.log");
check_file($out_file4, data_folder().$name."_out4.GSvar");

//varscan2 (tumor only) and COSMIC annotation
$out_file5 = output_folder().$name."_out5.GSvar";
check_exec("php ".src_folder()."/Tools/{$name}.php -in ".data_folder().$name."_in5.vcf -out $out_file5 -t_col GS150344_01 -db NGSD_TEST --log ".output_folder().$name."_out5.log");
check_file($out_file5, data_folder().$name."_out5.GSvar");

//input file contains VICC interpretation
$out_file6 = output_folder().$name."_out6.GSvar";
check_exec("php ".src_folder()."/Tools/{$name}.php -in ".data_folder().$name."_in6.vcf -out $out_file6 -t_col DX000001_01 -n_col DX000002_01 -db NGSD_TEST --log ".output_folder().$name."_out6.log");
check_file($out_file6, data_folder().$name."_out6.GSvar");

//umiVar2 (cfDNA)
$out_file7 = output_folder().$name."_out7.GSvar";
check_exec("php ".src_folder()."/Tools/{$name}.php -in ".data_folder().$name."_in7.vcf -out $out_file7 -t_col DXcf000001_01 -cfdna -db NGSD_TEST --log ".output_folder().$name."_out7.log");
check_file($out_file7, data_folder().$name."_out7.GSvar");
end_test();

?>
