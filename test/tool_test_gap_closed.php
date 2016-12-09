<?php

require_once("framework.php");

$name = "gap_closed";
$file = data_folder().$name;

start_test($name);

//output test
$out_file = output_folder().$name."_out1.tsv";
check_exec("php ".src_folder()."/Primer/".$name.".php -ab1 ".$file."_in1.ab1 ".$file."_in2.ab1 -gap ".$file."_in1.bed -out ".$out_file);
check_file($out_file, data_folder().$name."_out1.tsv");


$out_file = output_folder().$name."_out2.tsv";
check_exec("php ".src_folder()."/Primer/".$name.".php -ab1 ".$file."_in1.ab1 ".$file."_in2.ab1 -gap ".$file."_in1.bed -out ".$out_file." -hq -qual_thresh 23 -qual_window 2");
check_file($out_file, data_folder().$name."_out2.tsv");

end_test();

?>