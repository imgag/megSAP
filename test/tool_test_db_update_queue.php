<?php

require_once("framework.php");

$name = "db_update_queue";
$file = data_folder().$name;

start_test($name);

//init 
check_exec(get_path("ngs-bits")."NGSDInit -test -add ".data_folder()."/{$name}.sql");

//test - queued>started, cancel>canceled, running>finished, somatic is postponed
$out_file = output_folder().$name."_out1.txt";
check_exec("php ".src_folder()."/NGS/{$name}.php -db NGSD_TEST -debug > $out_file 2>&1");
check_file($out_file, data_folder().$name."_out1.txt");

//test 2 - running jobs finish, somatic is started
$out_file = output_folder().$name."_out2.txt";
check_exec("php ".src_folder()."/NGS/{$name}.php -db NGSD_TEST -debug > $out_file 2>&1");
check_file($out_file, data_folder().$name."_out2.txt");

end_test();

?>
