<?php

require_once("framework.php");
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/genomics.php");

$name = "rc_normalize";
start_test($name);

$in_file = data_folder().$name."_in.tsv";
$in_file_exon = data_folder().$name."_in_exon.tsv";
$prefix = output_folder().$name;

//test 1: only fpkm normalization
$out_file1 = output_folder().$name."_out1.tsv";
$log_file1 = output_folder().$name."_out1.log";
check_exec("php ".src_folder()."/NGS/{$name}.php -in '$in_file' -out '$out_file1' -method fpkm --log $log_file1");
check_file($out_file1, data_folder().$name."_out1.tsv");

//test 2: raw and cpm normalization
$out_file2 = output_folder().$name."_out2.tsv";
$log_file2 = output_folder().$name."_out2.log";
check_exec("php ".src_folder()."/NGS/{$name}.php -in '$in_file' -out '$out_file2' -method raw,cpm --log $log_file2");
check_file($out_file2, data_folder().$name."_out2.tsv");

//test 3: default parameters, raw, cpm, fpkm, tpm normalization
$out_file3 = output_folder().$name."_out3.tsv";
$log_file3 = output_folder().$name."_out3.log";
check_exec("php ".src_folder()."/NGS/{$name}.php -in '$in_file' -out '$out_file3' -method raw,cpm,fpkm,tpm --log $log_file3");
check_file($out_file3, data_folder().$name."_out3.tsv");

//test 4: exon-level
$out_file4 = output_folder().$name."_out4.tsv";
$out_file4_exon = output_folder().$name."_out4_exon.tsv";
$log_file4 = output_folder().$name."_out4.log";
check_exec("php ".src_folder()."/NGS/{$name}.php -in '$in_file' -in_exon '$in_file_exon' -out '$out_file4' -out_exon '$out_file4_exon' -method raw,cpm,fpkm,tpm --log $log_file3");
//gene-level: same result as test 3
check_file($out_file4, data_folder().$name."_out3.tsv");
check_file($out_file4_exon, data_folder().$name."_out4_exon.tsv");

end_test();

?>