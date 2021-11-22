<?php

require_once("framework.php");
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/genomics.php");

$name = "mapping_bwa";
start_test($name);

$in1_file = data_folder().$name."_in1.fastq.gz";
$in2_file = data_folder().$name."_in2.fastq.gz";

//test 1 (keep duplicates)
$out_file1 = output_folder().$name."_out1.bam";
$log_file1 = output_folder().$name."_out1.log";
check_exec("php ".src_folder()."/NGS/{$name}.php -in1 '$in1_file' -in2 '$in2_file' -out '$out_file1' --log $log_file1 -threads 1");
check_file($out_file1, data_folder().$name."_out1.bam");

//test 2 (remove dupliates)
$out_file2 = output_folder().$name."_out2.bam";
$log_file2 = output_folder().$name."_out2.log";
check_exec("php ".src_folder()."/NGS/{$name}.php -in1 '$in1_file' -in2 '$in2_file' -out '$out_file2' --log $log_file2 -dedup -threads 1");
check_file($out_file2, data_folder().$name."_out2.bam");

end_test();

?>