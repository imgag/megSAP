<?php

require_once("framework.php");

$name = "barcode_correction";
start_test($name);

$out_file1 = output_folder().$name."_out1.bam";
check_exec("python ".src_folder()."/NGS/barcode_correction.py --infile ".data_folder().$name."_in1.bam --outfile $out_file1 --BCerror 1 --minBQ 30");
check_file($out_file1, data_folder().$name."_out1.bam");

$out_file2 = output_folder().$name."_out2.bam";
check_exec("python ".src_folder()."/NGS/barcode_correction.py --infile ".data_folder().$name."_in2.bam --outfile $out_file2 --n --BCerror 1 --minBQ 30");
check_file($out_file2, data_folder().$name."_out2.bam");

end_test();

?>
