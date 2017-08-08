<?php

require_once("framework.php");
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/genomics.php");

$name = "rc_normalize";
start_test($name);

$in_file = data_folder().$name."_in.tsv";
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

//test 3: default parameters, raw, cpm, fpkm normalization
$out_file3 = output_folder().$name."_out3.tsv";
$log_file3 = output_folder().$name."_out3.log";
check_exec("php ".src_folder()."/NGS/{$name}.php -in '$in_file' -out '$out_file3' --log $log_file3");
check_file($out_file3, data_folder().$name."_out3.tsv");

end_test();

?>