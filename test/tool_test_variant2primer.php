<?php

require_once("framework.php");

$name = "variant2primer";
start_test($name);

//test
$in1 = data_folder().$name."_in1.txt";
$out1 = output_folder().$name."_out1.txt";
check_exec("php ".src_folder()."/Primer/".$name.".php -in $in1 -out $out1");

remove_lines_containing($out1, array("PRIMER_THERMODYNAMIC_PARAMETERS_PATH"));
check_file($out1, data_folder().$name."_out1.txt");

end_test();

?>
