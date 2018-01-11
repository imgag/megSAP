<?php

require_once("framework.php");

$name = "barcode_to_header";
start_test($name);

// standard strelka
$out_file1 = output_folder().$name."_out1.fastq.gz";
check_exec("python ".src_folder()."/NGS/barcode_to_header.py -i ".data_folder().$name."_in1.fastq.gz -o $out_file1 -bc1 ".data_folder().$name."_in2.fastq.gz -bc2 ".data_folder().$name."_in3.fastq.gz");
check_file($out_file1, data_folder().$name."_out1.fastq.gz");

end_test();

?>
