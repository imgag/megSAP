<?php

require_once("framework.php");

$name = "converter_tsv2bed";
start_test($name);

$in_file1 = data_folder().$name."_in.tsv";
$out_file = output_folder().$name."_out_test.bed";
$reference_file = data_folder().$name."_out.bed";

check_exec("php ".src_folder()."/Tools/converter_tsv2bed.php -tsv '$in_file1' -out '$out_file'");
check_file($out_file, $reference_file);

end_test();

?>