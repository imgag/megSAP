<?php

require_once("framework.php");
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/genomics.php");

$name = "read_counting_featureCounts";
start_test($name);

$in = data_folder().$name."_in.bam";
$prefix = output_folder().$name;

$out_file = output_folder().$name."_out.tsv";
$log_file = output_folder().$name."_out.log";

check_exec("php ".src_folder()."/NGS/{$name}.php -in '$in' -prefix '$prefix' -paired --log $log_file");
// Change to canonical output filename
rename("${prefix}_counts.tsv", $out_file);
check_file($out_file, data_folder().$name."_out.tsv");

end_test();

?>