<?php

require_once("framework.php");

$name = "db_add_mids";
start_test($name);

//init
check_exec(get_path("ngs-bits")."NGSDInit -test");

//test
$log_file = output_folder()."/{$name}_out1.log";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder()."/{$name}.tsv -db NGSD_TEST --log $log_file");

end_test();

?>