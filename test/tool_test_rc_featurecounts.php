<?php

require_once("framework.php");
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/genomics.php");

$name = "rc_featurecounts";
start_test($name);


$out = output_folder().$name."_out1.txt";
check_exec("php ".src_folder()."/Tools/{$name}.php -in ".data_folder().$name."_in.bam -out $out -library_type unstranded --log ".output_folder().$name."_out1.log");
remove_lines_containing($out, array("Program:featureCounts", "rc_featurecounts_in.bam"));
exec2("gzip $out"); //gzip output - test data would be very large otherwise
check_file($out.".gz", data_folder().$name."_out1.tsv.gz", false);

end_test();

?>