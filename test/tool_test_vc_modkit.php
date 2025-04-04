<?php

require_once("framework.php");

$name = "vc_modkit";
start_test($name);

########################## run modkit with summary ##################
check_exec("php ".src_folder()."/Tools/{$name}.php -bam ".data_folder().$name."_in1.bam -bed ".output_folder().$name."_out1.bed.gz -summary ".output_folder().$name."_summary_out1.txt --log ".output_folder().$name."_out1.log");
check_file(output_folder().$name."_out1.bed.gz", data_folder().$name."_out1.bed.gz");

end_test();

?>


