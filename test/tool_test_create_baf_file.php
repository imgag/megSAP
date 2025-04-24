<?php

require_once("framework.php");

$name = "create_baf_file";
start_test($name);

//no variants left to score after filtering > input is equal to output
$out1 = output_folder().$name."_out1.tsv";
check_exec("php ".src_folder()."/Auxilary/{$name}.php -vcf ".data_folder().$name."_in1.vcf.gz -s_col DNA2501927A1_01 -out_file $out1 --log ".output_folder().$name."_out1.log");
check_file($out1, data_folder().$name."_out1.tsv", true);

end_test();

?>
