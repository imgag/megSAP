<?php

require_once("framework.php");
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/genomics.php");

$name = "rc_annotate";
start_test($name);

// Annotate Ensembl gene identifiers with HGNC gene identifier/gene names
$in_file1 = data_folder().$name."_in1.tsv";
$out_file1 = output_folder().$name."_out1.tsv";
check_exec("php ".src_folder()."/NGS/{$name}.php -in '$in_file1' -out '$out_file1'");
check_file($out_file1, data_folder().$name."_out1.tsv");

// Annotate Ensembl gene identifiers with HGNC gene identifier/gene names
$in_file2 = data_folder().$name."_in2.tsv";
$out_file2 = output_folder().$name."_out2.tsv";
check_exec("php ".src_folder()."/NGS/{$name}.php -in '$in_file2' -out '$out_file2'");
check_file($out_file2, data_folder().$name."_out2.tsv");

end_test();

?>