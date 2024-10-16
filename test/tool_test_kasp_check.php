<?php

require_once("framework.php");

$name = "kasp_check";
start_test($name);

if (db_is_enabled("NGSD"))
{
	//test 1 - LightCycler - SNP set 1
	$in1 = data_folder().$name."_in1.txt";
	$out1 = output_folder().$name."_out1.tsv";

	$stdout = check_exec("php ".src_folder()."/Tools/".$name.".php -in $in1 -snps set1 -format LightCycler -out $out1 -user unknown");
	unlink(data_folder().$name."_in1_converted.tsv");

	check_file($out1, data_folder().$name."_out1.tsv");


	//test 2 - LightCycler - SNP set 2
	$in2 = data_folder().$name."_in2.txt";
	$out2 = output_folder().$name."_out2.tsv";

	$stdout = check_exec("php ".src_folder()."/Tools/".$name.".php -in $in2 -snps set2 -format LightCycler -out $out2 -user unknown");
	unlink(data_folder().$name."_in2_converted.tsv");

	check_file($out2, data_folder().$name."_out2.tsv");

	
	//test 3 - StepOnePlus - SNP set 2 - only conversion is checked
	$in3 = data_folder().$name."_in3.txt";
	$out3 = output_folder().$name."_out3.tsv";
	$converted3 = output_folder().$name."_converted3.tsv";

	$stdout = check_exec("php ".src_folder()."/Tools/".$name.".php -in $in3 -snps set2 -format StepOnePlus -out $out3 -user unknown");
	rename(data_folder().$name."_in3_converted.tsv", $converted3); 
	check_file($converted3, data_folder().$name."_converted3.tsv");
}

end_test();

?>