<?php

require_once("framework.php");
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/genomics.php");

$name = "mapping_minimap";
start_test($name);

$in1_file = data_folder().$name."_in1.fastq.gz";
$in2_file = data_folder().$name."_out1.bam";
$in1_system =  data_folder().$name."_in1_system.ini";

//test 1 (keep duplicates)
$out_file1 = output_folder().$name."_out1.bam";
check_exec("php ".src_folder()."/NGS/{$name}.php -in_fastq {$in1_file} -system {$in1_system} -out {$out_file1} --log ".output_folder().$name."_out1.log -threads 1");
check_file($out_file1, data_folder().$name."_out1.bam");

//bam input
$out_file2 = output_folder().$name."_out2.bam";
check_exec("php ".src_folder()."/NGS/{$name}.php -in_bam {$in2_file} -system {$in1_system} -out {$out_file2} --log ".output_folder().$name."_out2.log -threads 1");
check_file($out_file2, data_folder().$name."_out1.bam");


end_test();

?>
