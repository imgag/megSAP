<?php

require_once("framework.php");

$name = "converter_keg2tsv";
start_test($name);

$in1 = data_folder().$name."_in1.keg";
$out1 = output_folder().$name."_out1.tsv";
check_exec("php ".src_folder()."/Tools/".$name.".php -in $in1 -out $out1");

check_file($out1, data_folder().$name."_out1.tsv");

end_test();

?>