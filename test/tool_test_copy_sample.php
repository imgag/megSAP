<?php

require_once("framework.php");

$name = "copy_sample";
$file = data_folder().$name;

start_test($name);

if (db_is_enabled("NGSD"))
{
	mkdir("Unaligned");

	$out_file = output_folder().$name."_out1Makefile";
	check_exec("php ".src_folder()."/NGS/".$name.".php -samplesheet ".$file."_in1.csv -out ".$out_file);
	check_file($out_file, data_folder().$name."_out1Makefile");

	$out_file = output_folder().$name."_out2Makefile";
	check_exec("php ".src_folder()."/NGS/".$name.".php -high_priority -overwrite -samplesheet ".$file."_in2.csv -out ".$out_file);
	check_file($out_file, data_folder().$name."_out2Makefile");

	rmdir("Unaligned");
}

end_test();

?>
