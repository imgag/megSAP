<?php

require_once("framework.php");

$name = "region2exon";
start_test($name);

//test 1
$out1 = output_folder().$name."_out1.txt";
check_exec("php ".src_folder()."/Primer/".$name.".php -reg chr1:151137550-151138268 -out $out1");
check_file($out1, data_folder().$name."_out1.txt");

//test 2
$out2 = output_folder().$name."_out2.txt";
check_exec("php ".src_folder()."/Primer/".$name.".php -reg chr17:41244361-41244377 -out $out2");
check_file($out2, data_folder().$name."_out2.txt");

//test 3
$out3 = output_folder().$name."_out3.txt";
check_exec("php ".src_folder()."/Primer/".$name.".php -reg chr1:151139409-151139498 -out $out3");
check_file($out3, data_folder().$name."_out3.txt");

//test 4
$out4 = output_folder().$name."_out4.txt";
check_exec("php ".src_folder()."/Primer/".$name.".php -reg chr17:40044361-40044377 -out $out4");
check_file($out4, data_folder().$name."_out4.txt");

end_test();

?>