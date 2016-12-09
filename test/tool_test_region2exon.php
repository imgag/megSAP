<?php

require_once("framework.php");

$name = "region2exon";
start_test($name);

//test
$in1 = data_folder().$name."_in.txt";
$out1 = output_folder().$name."_out.txt";
check_exec("php ".src_folder()."/Primer/".$name.".php -in $in1 -out $out1");
print "php ".src_folder()."/Primer/".$name.".php -in $in1 -out $out1";
check_file($out1, data_folder().$name."_out.txt");

end_test();

?>