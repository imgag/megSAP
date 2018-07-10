<?php

require_once("framework.php");

$name = "filter_bam";

start_test($name."_1");
$out_file1 = output_folder().$name."_out1.bam";
check_exec(get_path("ngs-bits")."BamFilter -in ".data_folder().$name."_in1.bam -out $out_file1 -minMQ 30 -maxMM 4 -maxGap 2");
check_file($out_file1, data_folder().$name."_out1.bam");
end_test();

start_test($name."_2");
$out_file2 = output_folder().$name."_out2.bam";
check_exec(get_path("ngs-bits")."BamFilter -in ".data_folder().$name."_in2.bam -out $out_file2 -minMQ 0 -maxMM 20 -maxGap 20 -minDup 2");
check_file($out_file2, data_folder().$name."_out2.bam");
end_test();

?>
