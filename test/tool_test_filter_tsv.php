<?php

require_once("framework.php");

$name = "filter_tsv";

start_test($name);

$out_file1 = output_folder().$name."_out1.tsv";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in1.tsv -out $out_file1 -type coding");
check_file($out_file1, data_folder().$name."_out1.tsv");

$out_file2 = output_folder().$name."_out2.tsv";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in1.tsv -out $out_file2 -type somatic");
check_file($out_file2, data_folder().$name."_out2.tsv");

end_test();
?>
