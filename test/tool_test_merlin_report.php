<?php

require_once("framework.php");

$name = "merlin_report";
$file = data_folder().$name;

start_test($name);

//output test
$out_file = output_folder().$name."_out1.bed";
check_exec("php ".src_folder()."/Chips/".$name.".php -in ".$file."_in1.tbl -out ".$out_file);
check_file($out_file, data_folder().$name."_out1.bed");

end_test();

?>