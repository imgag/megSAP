<?php

require_once("framework.php");

$name = "an_aidiva";
start_test($name);

//standard
$out_file1 = output_folder().$name."_out1.vcf";
check_exec("php ".src_folder()."/Tools/{$name}.php -in ".data_folder().$name."_in1.vcf -out $out_file1 --log ".output_folder().$name."_out1.log");
check_file($out_file1, data_folder().$name."_out1.vcf", true);

?>
