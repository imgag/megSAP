<?php

require_once("framework.php");

$name = "barcode_correction";
start_test($name);

$out_file1 = output_folder().$name."_out1.bam";
check_exec(get_path("python3")." ".get_path("umiVar2")."/barcode_correction.py --infile ".data_folder().$name."_in1.bam --outfile $out_file1 --barcode_error 1 --minBQ 30");
check_file($out_file1, data_folder().$name."_out1.bam");

$out_file2 = output_folder().$name."_out2.bam";
check_exec(get_path("python3")." ".get_path("umiVar2")."/barcode_correction.py --infile ".data_folder().$name."_in2.bam --outfile $out_file2 --n --barcode_error 1 --minBQ 30");
check_exec("samtools sort -o $out_file2 $out_file2");
check_file($out_file2, data_folder().$name."_out2.bam");

end_test();

?>
