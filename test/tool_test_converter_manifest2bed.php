<?php

require_once("framework.php");

$name = "converter_manifest2bed";
start_test($name);

//
$in_file = data_folder().$name."_in.txt";
$out_file = output_folder().$name."_out.bed";
check_exec("php ".src_folder()."/Tools/converter_manifest2bed.php -m $in_file -o $out_file");
check_file($out_file, data_folder().$name."_out.bed");


end_test();
?>
