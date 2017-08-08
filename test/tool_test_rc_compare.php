<?php

require_once("framework.php");
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/genomics.php");

$name = "rc_compare";
start_test($name);

$in_file1 = data_folder().$name."_in1.tsv";
$in_file2 = data_folder().$name."_in2.tsv";

$out_file1 = output_folder().$name."_out1.tsv";
check_exec("php ".src_folder()."/NGS/{$name}.php -in1 '$in_file1' -in2 '$in_file2' -out '$out_file1'");
check_file($out_file1, data_folder().$name."_out1.tsv");

$out_file2 = output_folder().$name."_out2.tsv";
check_exec("php ".src_folder()."/NGS/{$name}.php -in1 '$in_file1' -in2 '$in_file2' -out '$out_file2' -column cpm -method log1p_fc");
check_file($out_file2, data_folder().$name."_out2.tsv");

$out_file3 = output_folder().$name."_out3.tsv";
check_exec("php ".src_folder()."/NGS/{$name}.php -in1 '$in_file1' -in2 '$in_file2' -out '$out_file3' -column cpm -method log_fc");
check_file($out_file3, data_folder().$name."_out3.tsv");

end_test();

?>