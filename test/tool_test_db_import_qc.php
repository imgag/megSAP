<?php

require_once("framework.php");

$name = "db_import_qc";
start_test($name);

//init
check_exec(get_path("ngs-bits")."NGSDInit -test -add ".data_folder()."/{$name}.sql");

//import 3 samples
$log1_file = output_folder().$name."_out1.log";
check_exec("php ".src_folder()."/NGS/db_import_qc.php -id GS120676_01 -files ".data_db_folder()."GS120676_01*.qcML -db NGSD_TEST --log $log1_file");

$log2_file = output_folder().$name."_out2.log";
check_exec("php ".src_folder()."/NGS/db_import_qc.php -id GS120677_01 -files ".data_db_folder()."GS120677_01*.qcML -db NGSD_TEST --log $log2_file");

$log3_file = output_folder().$name."_out3.log";
check_exec("php ".src_folder()."/NGS/db_import_qc.php -id GS120678_01 -files ".data_db_folder()."GS120678_01*.qcML -db NGSD_TEST --log $log3_file");

//import third sample again (with -force flag)
$log4_file = output_folder().$name."_out4.log";
check_exec("php ".src_folder()."/NGS/db_import_qc.php -id GS120678_01 -files ".data_db_folder()."GS120678_01*.qcML -db NGSD_TEST --log $log4_file -force");


end_test();

?>