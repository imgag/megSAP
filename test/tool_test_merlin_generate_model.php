<?php

require_once("framework.php");

$name = "merlin_generate_model";
$file = data_folder().$name;

start_test($name);

//output test
$out_file = output_folder().$name."_out1.model";
check_exec("php ".src_folder()."/Chips/".$name.".php -type dominant_hp -dat ".$file."_in1.dat -ped ".$file."_in1.ped -out ".$out_file);
check_file($out_file, data_folder().$name."_out1.model");

end_test();

?>