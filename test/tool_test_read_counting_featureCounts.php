<?php

require_once("framework.php");
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/genomics.php");

$name = "read_counting_featureCounts";
start_test($name);

$in = data_folder().$name."_in.bam";
$prefix = output_folder().$name;
$out = $prefix."_counts.tsv";

$log_file = output_folder().$name."_out.log";

check_exec("php ".src_folder()."/NGS/{$name}.php -in '$in' -prefix '$prefix' -paired --log $log_file");
remove_lines_containing($out, array("Program:featureCounts", "read_counting_featureCounts_in.bam"));
check_file($out, data_folder().$name."_out.tsv", false);

end_test();

?>