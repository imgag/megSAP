<?php

require_once("framework.php");
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/genomics.php");

$name = "annotate_count_file";
start_test($name);

$in_file = data_folder().$name."_in.tsv";

// Annotate RefSeq identifier with HGNC gene identifier
$out_file = output_folder().$name."_out.tsv";
$log_file = output_folder().$name."_out.log";
check_exec("php ".src_folder()."/NGS/{$name}.php -in '$in_file' -out '$out_file'");
check_file($out_file, data_folder().$name."_out.tsv");

end_test();

?>