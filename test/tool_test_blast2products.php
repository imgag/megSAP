<?php

require_once("framework.php");

$name = "blast2products";
start_test($name);

//test
$in1 = data_folder().$name."_in1.blast";
$in2 = data_folder().$name."_in2.blast";
$out1 = output_folder().$name."_out.txt";
check_exec("php ".src_folder()."/Primer/".$name.".php -in $in1 -in2 $in2 -out $out1");
print "php ".src_folder()."/Primer/".$name.".php -in $in1 -in2 $in2 -out $out1";
check_file($out1, data_folder().$name."_out.txt");

end_test();

?>