<?php

require_once("framework.php");

$name = "kasp_check";
start_test($name);

if (production_ngsd_enabled())
{
	//test 1 - SNP set 1
	$in1 = data_folder().$name."_in1.txt";
	$out1 = output_folder().$name."_out1.tsv";

	$stdout = check_exec("php ".src_folder()."/Tools/".$name.".php -in $in1 -snps set1 -out $out1");
	unlink(data_folder().$name."_in1_converted.tsv");

	check_file($out1, data_folder().$name."_out1.tsv");


	//test 2 - SNP set 2
	$in2 = data_folder().$name."_in2.txt";
	$out2 = output_folder().$name."_out2.tsv";

	$stdout = check_exec("php ".src_folder()."/Tools/".$name.".php -in $in2 -snps set2 -out $out2");

	check_file($out2, data_folder().$name."_out2.tsv");
}

end_test();

?>