<?php

require_once("framework.php");

$name = "converter_bgi2tsv";
start_test($name);

$in_file1 = data_folder().$name."_in.gff";
$out_file = output_folder().$name."_out_test.tsv";
$reference_file = data_folder().$name."_out.tsv";

check_exec("php ".src_folder()."/Tools/converter_bgi2tsv.php -gff '$in_file1' -out '$out_file'");
check_file($out_file, $reference_file);

end_test();

?>