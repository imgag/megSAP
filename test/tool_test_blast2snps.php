<?php

require_once("framework.php");

$name = "blast2snps";
start_test($name);

//test
$out_file = output_folder().$name."_out1.txt";
check_exec("php ".src_folder()."/Primer/".$name.".php -in ".data_folder().$name."_in1.txt -blast ".data_folder().$name."_blast_in1.txt -out $out_file -db_list dbSNP,ESP6500,ExAC");
check_file($out_file, data_folder().$name."_out1.txt");

end_test();

?>