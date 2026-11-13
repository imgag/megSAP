<?php

require_once("framework.php");

$name = "vc_modkit";
start_test($name);

// cram test:
check_exec("php ".src_folder()."/Tools/{$name}.php -threads 4 -region chr20 -bam ".data_folder().$name."_in2.cram -bed ".output_folder().$name."_out2.bed.gz -summary ".output_folder().$name."_summary_out2.txt --log ".output_folder().$name."_out2.log");
check_file(output_folder().$name."_out2.bed.gz", data_folder().$name."_out2.bed.gz");
check_file(output_folder().$name."_out2_hp1.bed.gz", data_folder().$name."_out2_hp1.bed.gz");
check_file(output_folder().$name."_out2_hp2.bed.gz", data_folder().$name."_out2_hp2.bed.gz");
check_file(output_folder().$name."_summary_out2.txt", data_folder().$name."_summary_out2.txt");

// bam test:
check_exec("php ".src_folder()."/Tools/{$name}.php -threads 4 -region chr20 -bam ".data_folder().$name."_in2.cram -bed ".output_folder().$name."_out3.bed.gz -summary ".output_folder().$name."_summary_out3.txt --log ".output_folder().$name."_out3.log");
check_file(output_folder().$name."_out3.bed.gz", data_folder().$name."_out2.bed.gz");
check_file(output_folder().$name."_out3_hp1.bed.gz", data_folder().$name."_out2_hp1.bed.gz");
check_file(output_folder().$name."_out3_hp2.bed.gz", data_folder().$name."_out2_hp2.bed.gz");
check_file(output_folder().$name."_summary_out3.txt", data_folder().$name."_summary_out2.txt");


end_test();

?>


