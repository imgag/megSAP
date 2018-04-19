<?php

require_once("framework.php");

$name = "db_queue_analysis";
$file = data_folder().$name;

start_test($name);

//init
check_exec(get_path("ngs-bits")."NGSDInit -test -add ".data_folder()."/{$name}.sql");

//test single
$log_file = output_folder().$name."_out1.log";
check_exec("php ".src_folder()."/NGS/{$name}.php -user ahklauo1 -db NGSD_TEST -type 'single sample' -samples DX181277_01 -args '-bla' --log $log_file");

//test multi
$log_file = output_folder().$name."_out2.log";
check_exec("php ".src_folder()."/NGS/{$name}.php -user ahklauo1 -db NGSD_TEST -type 'multi sample' -samples DX181277_01 DX181278_01 DX181279_01 -info affected control control -args '-bli -bla -bluff' --log $log_file");

//test trio
$log_file = output_folder().$name."_out3.log";
check_exec("php ".src_folder()."/NGS/{$name}.php -user ahklauo1 -db NGSD_TEST -type 'trio' -samples DX181277_01 DX181278_01 DX181279_01 -info child father mother -high_priority --log $log_file");

//test somatic
$log_file = output_folder().$name."_out4.log";
check_exec("php ".src_folder()."/NGS/{$name}.php -user ahklauo1 -db NGSD_TEST -type 'somatic' -samples DX181277_01 DX181278_01 -info tumor normal -high_priority --log $log_file");

end_test();

?>
