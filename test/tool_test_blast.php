<?php

require_once("framework.php");

$name = "blast";
start_test($name);

//test
$out_file = output_folder().$name."_out1.txt";
check_exec("php ".src_folder()."/Primer/".$name.".php -in ".data_folder().$name."_in1.fasta -out $out_file");
check_file($out_file, data_folder().$name."_out1.txt");

end_test();

?>