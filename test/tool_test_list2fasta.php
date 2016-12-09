<?php

require_once("framework.php");

$name = "list2fasta";
start_test($name);

//test
$out_file = output_folder().$name."_out1.fasta";
check_exec("php ".src_folder()."/Primer/".$name.".php -in ".data_folder().$name."_in1.txt -out $out_file");
check_file($out_file, data_folder().$name."_out1.fasta");

end_test();

?>