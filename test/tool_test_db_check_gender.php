<?php

require_once("framework.php");

$name = "db_check_gender";
start_test($name);

$bam = get_path("test_data_folder")."/GS140127_01.bam";
if (file_exists($bam))	
{

	//test with gender from CLI
	$log_file = output_folder().$name."_out1.log";
	check_exec("php ".src_folder()."/NGS/db_check_gender.php -in $bam -pid GS140127_01 -gender male --log $log_file");

	//test with gender from NGSD
	init_ngsd($name);
	$log_file = output_folder().$name."_out2.log";
	check_exec("php ".src_folder()."/NGS/db_check_gender.php -in $bam -pid GS140127_01 -db NGSD_TEST --log $log_file");
}

end_test();

?>