<?php

require_once("framework.php");

$name = "filter_confidence";

start_test($name);

//threshold 0.1 => no filtering
$outfile = output_folder().$name."_out1.txt";
check_exec("php ".src_folder()."/Chips/".$name.".php -in ".data_folder().$name."_in.txt -out $outfile -thres 0.1");
check_file($outfile, data_folder().$name."_out.txt");

//threshold 0.02 => filtering
$outfile = output_folder().$name."_out2.txt";
check_exec("php ".src_folder()."/Chips/".$name.".php -in ".data_folder().$name."_in.txt -out $outfile -thres 0.02");
check_file($outfile, data_folder().$name."_out2.txt");

end_test();

?>