<?php

require_once("framework.php");

$name = "primer2SNP";
start_test($name);

//test
$out_file = output_folder().$name."_out.txt";
check_exec("php ".src_folder()."/Primer/".$name.".php -in ".data_folder().$name."_in.bed -out $out_file");
check_file($out_file, data_folder().$name."_out.txt");

//test
$out_file = output_folder().$name."_out2.txt";
check_exec("php ".src_folder()."/Primer/".$name.".php -in ".data_folder().$name."_in2.bed -out $out_file -min_freq 0.0001");
check_file($out_file, data_folder().$name."_out2.txt");


end_test();

?>