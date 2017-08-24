<?php

require_once("framework.php");
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/genomics.php");

$name = "mapping_star";
start_test($name);

$in1_file = data_folder().$name."_in1.fastq.gz";
$in2_file = data_folder().$name."_in2.fastq.gz";

//test 1 (keep duplicates)
$out_file1 = output_folder().$name."_out1.bam";
check_exec("php ".src_folder()."/NGS/{$name}.php -in1 $in1_file -in2 $in2_file -out $out_file1 --log ".output_folder().$name."_out1.log -threads 1 -skip_dedup");
check_file($out_file1, data_folder().$name."_out1.bam");

end_test();

?>
