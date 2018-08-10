<?php

require_once("framework.php");

$name = "compare_array";
start_test($name);

//test1
$out_file1 = output_folder().$name."_out1.txt";
$out_file1 = output_folder().$name."_out1.txt";
check_exec("php ".src_folder()."/Tools/{$name}.php -a ".data_folder().$name."_cytoscan750K.txt -n ".data_folder().$name."_ngs.GSvar > $out_file1 2>&1");
check_file($out_file1, data_folder().$name."_out1.txt", true);

end_test();

?>
