<?php

require_once("framework.php");
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/genomics.php");

$name = "mapping_blatlike";
start_test($name);

//test 1
$out_file1 = output_folder().$name."_out1.tsv";
$log_file1 = output_folder().$name."_out1.log";
check_exec("php ".src_folder()."/Tools/{$name}.php -seq CTGGGTGACAGAGTGAGACCCTGTCTCTTTAAAACAATAGAAGTAGGCCAGGCACAGTGGCTCATGCCTGTAATCCAAGCACTTTGGGAGGTCGAGGCGGGTGGATCACGAGGTCAGGAGCTCCAGACCATCCTGGCTAACACGGTAAAACCCCGTCT -secondary 4 -out $out_file1 --log $log_file1");
check_file($out_file1, data_folder().$name."_out1.tsv");

end_test();

?>