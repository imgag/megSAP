<?php

require_once("framework.php");

$name = "db_check_gender";
start_test($name);

//init
check_exec(get_path("ngs-bits")."NGSDInit -test -add ".data_folder()."/{$name}.sql");

//test
$log_file = output_folder().$name."_out1.log";
check_exec("php ".src_folder()."/NGS/db_check_gender.php -in ".get_path("data_folder")."/test_data/GS120159.bam -pid GS120159_01 -db NGSD_TEST --log $log_file");

//check that 'male' was detected
$is_male = false;
foreach(file($log_file) as $line)
{
	if (contains($line, "gender: male")) $is_male = true;
}
check($is_male, true);

end_test();

?>