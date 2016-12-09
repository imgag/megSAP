<?php

require_once("framework.php");

$name = "plink_regions";
start_test($name);

$out_file = output_folder().$name."_out.txt";
check_exec("php ".src_folder()."/Chips/".$name.".php -in ".data_folder().$name."_in.txt -gap 1000000 -out $out_file");
check_file($out_file, data_folder().$name."_out.txt");

$out_file = output_folder().$name."_out2.txt";
check_exec("php ".src_folder()."/Chips/".$name.".php -in ".data_folder().$name."_in2.txt -gap 1000000 -out $out_file");
check_file($out_file, data_folder().$name."_out2.txt");

end_test();

?>