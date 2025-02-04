<?php

require_once("framework.php");

$name = "filter_kraken2";
start_test($name);

$in1 = data_folder().$name."_in1.fastq.gz";
$in2 = data_folder().$name."_in2.fastq.gz";

$out_files = output_folder().$name."_out#.fastq";
check_exec("php ".src_folder()."/Tools/{$name}.php -output - -paired -gzip_compressed -unclassified_out $out_files -in $in1 $in2 --log ".output_folder().$name.".log");
check_file(output_folder().$name."_out_1.fastq", data_folder().$name."_out1.fastq");
check_file(output_folder().$name."_out_2.fastq", data_folder().$name."_out2.fastq");

end_test();
?>
