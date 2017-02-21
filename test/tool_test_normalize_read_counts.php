<?php

require_once("framework.php");
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/genomics.php");

$name = "normalize_read_counts";
start_test($name);

$in_file = data_folder().$name."_in.tsv";
$prefix = output_folder().$name;

//test 1: rpkm normalization
$out_file1 = output_folder().$name."_out1.tsv";
$log_file1 = output_folder().$name."_out1.log";
check_exec("php ".src_folder()."/NGS/{$name}.php -in '$in_file' -out '$out_file1' -method rpkm -header --log $log_file1");
check_file($out_file1, data_folder().$name."_out1.tsv");

//test 2: cpm normalization
$out_file2 = output_folder().$name."_out2.tsv";
$log_file2 = output_folder().$name."_out2.log";
check_exec("php ".src_folder()."/NGS/{$name}.php -in '$in_file' -out '$out_file2' -method cpm -header --log $log_file2");
check_file($out_file2, data_folder().$name."_out2.tsv");

end_test();

?>