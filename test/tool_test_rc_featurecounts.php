<?php

require_once("framework.php");
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/genomics.php");

$name = "rc_featurecounts";
start_test($name);


$out = output_folder().$name."_out1.tsv";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in.bam -out $out -library_type unstranded --log ".output_folder().$name."_out1.log");
remove_lines_containing($out, array("Program:featureCounts", "rc_featurecounts_in.bam"));
check_file($out, data_folder().$name."_out1.tsv", false);

end_test();

?>