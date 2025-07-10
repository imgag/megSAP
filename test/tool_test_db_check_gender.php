<?php

require_once("framework.php");

$name = "db_check_gender";
start_test($name);

$bam = get_path("test_data_folder")."/GS140127_01.bam";
if (file_exists($bam))	
{

	//test with gender from CLI
	$log_file = output_folder().$name."_out1.log";
	check_exec("php ".src_folder()."/Tools/db_check_gender.php -in $bam -pid GS140127_01 -gender male --log $log_file");

	//test with gender from NGSD
	init_ngsd($name);
	$log_file = output_folder().$name."_out2.log";
	check_exec("php ".src_folder()."/Tools/db_check_gender.php -in $bam -pid GS140127_01 -db NGSD_TEST --log $log_file");
	
	//test also checking SRY coverage
	$log_file = output_folder().$name."_out3.log";
	check_exec("php ".src_folder()."/Tools/db_check_gender.php -in $bam -pid GS140127_01 -gender male -check_sry_cov --log $log_file");

	//check failing
	$log_file = output_folder().$name."_out4.log";
	check_exec("php ".src_folder()."/Tools/db_check_gender.php -in $bam -pid GS140127_01 -gender female --log $log_file", false);
	
	$error_found = false;
	$lines = file($log_file);
	foreach($lines as $line)
	{
		if (strpos($line, "Expected gender 'female', but determined gender 'male' using method hetx on BAM file") !== false)
		{
			$error_found = true;
		}
	}

	check($error_found, true);
	


}

end_test();

?>