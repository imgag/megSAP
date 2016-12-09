<?php

require_once("framework.php");

$name = "bed_annotation";
start_test($name);

$out_file = output_folder().$name."_out1.tsv";
check_exec("php ".src_folder()."/Chips/".$name.".php -in ".data_folder().$name."_in1.bed -out $out_file");
check_file($out_file, data_folder().$name."_out1.tsv");

$out_file = output_folder().$name."_out2.tsv";
check_exec("php ".src_folder()."/Chips/".$name.".php -in ".data_folder().$name."_in2.bed -out $out_file");
check_file($out_file, data_folder().$name."_out2.tsv");

end_test();

?>