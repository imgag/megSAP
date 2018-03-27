<?php

require_once("framework.php");

$name = "filter_bam";

start_test($name."_1");
$out_file1 = output_folder().$name."_out1.bam";
check_exec("python ".src_folder()."/NGS/filter_bam.py --infile ".data_folder().$name."_in1.bam --outfile $out_file1 --minMQ 30 --maxMM 4 --maxGAP 2");
check_file($out_file1, data_folder().$name."_out1.bam");
end_test();

start_test($name."_2");
$out_file2 = output_folder().$name."_out2.bam";
check_exec("python ".src_folder()."/NGS/filter_bam.py --infile ".data_folder().$name."_in2.bam --outfile $out_file2 --minMQ 0 --maxMM 20 --maxGAP 20 --minDUP 2");
check_file($out_file2, data_folder().$name."_out2.bam");
end_test();

?>
