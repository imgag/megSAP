<?php

require_once("framework.php");

$name = "db_check_gender";
start_test($name);

//test with gender from CLI
$log_file = output_folder().$name."_out1.log";
check_exec("php ".src_folder()."/NGS/db_check_gender.php -in ".get_path("test_data_folder")."/GS140127_01.bam -pid GS140127_01 -gender male --log $log_file");

//test with gender from NGSD
check_exec(get_path("ngs-bits")."NGSDInit -test -add ".data_folder()."/{$name}.sql");
$log_file = output_folder().$name."_out2.log";
check_exec("php ".src_folder()."/NGS/db_check_gender.php -in ".get_path("test_data_folder")."/GS140127_01.bam -pid GS140127_01 -db NGSD_TEST --log $log_file");

end_test();

?>