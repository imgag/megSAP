<?php

require_once("framework.php");
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/genomics.php");

$name = "mapping_star_htseq";
start_test($name);

$in1_file = data_folder().$name."_in1.fastq.gz";
$in2_file = data_folder().$name."_in2.fastq.gz";
$ini_file = data_folder().$name.".ini";

//test 1 (keep duplicates)
$prefix = output_folder().$name;
$out_file1 = output_folder().$name."_out1.bam";
$log_file1 = output_folder().$name."_out1.log";
check_exec("php ".src_folder()."/NGS/{$name}.php -in1 '$in1_file' -in2 '$in2_file' -prefix '$prefix' -system $ini_file --log $log_file1 -p 1");
// Change to canonical output filename
rename("${prefix}Aligned.sortedByCoord.out.bam", $out_file1);
check_file($out_file1, data_folder().$name."_out1.bam");

end_test();

?>