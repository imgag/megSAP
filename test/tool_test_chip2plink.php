<?php

require_once("framework.php");

$name = "chip2plink";
start_test($name);

// illumina 6k chip
$map = output_folder().$name."_out1_map.txt";
$ped = output_folder().$name."_out1_ped.txt";
$fam = output_folder().$name."_out1_fam.txt";

check_exec("php ".src_folder()."/Chips/".$name.".php -chip ".data_folder().$name."_in1.txt -meta ".data_folder().$name."_in1_meta.txt -type illumina_6k -family TEST -out_map $map -out_ped $ped -out_fam $fam");
check_file($fam, data_folder().$name."_out1_fam.txt");
check_file($ped, data_folder().$name."_out1_ped.txt");
check_file($map, data_folder().$name."_out1_map.txt");

// affymetrix 6.0 chip without confidence values
$map = output_folder().$name."_out1_map.txt";
$ped = output_folder().$name."_out1_ped.txt";
$fam = output_folder().$name."_out1_fam.txt";

check_exec("php ".src_folder()."/Chips/".$name.".php -chip ".data_folder().$name."_in2.txt -meta ".data_folder().$name."_in2_meta.txt -type affymetrix_6.0 -family TEST -out_map $map -out_ped $ped -out_fam $fam");
check_file($fam, data_folder().$name."_out2_fam.txt");
check_file($ped, data_folder().$name."_out2_ped.txt");
check_file($map, data_folder().$name."_out2_map.txt");

// affymetrix 6.0 chip with confidence values
$map = output_folder().$name."_out1_map.txt";
$ped = output_folder().$name."_out1_ped.txt";
$fam = output_folder().$name."_out1_fam.txt";

check_exec("php ".src_folder()."/Chips/".$name.".php -chip ".data_folder().$name."_in3.txt -meta ".data_folder().$name."_in2_meta.txt -type affymetrix_6.0 -family TEST -out_map $map -out_ped $ped -out_fam $fam");
check_file($fam, data_folder().$name."_out2_fam.txt");
check_file($ped, data_folder().$name."_out2_ped.txt");
check_file($map, data_folder().$name."_out2_map.txt");

// illumina cytoscan_hd chip
$map = output_folder().$name."_out4_map.txt";
$ped = output_folder().$name."_out4_ped.txt";
$fam = output_folder().$name."_out4_fam.txt";

check_exec("php ".src_folder()."/Chips/".$name.".php -chip ".data_folder().$name."_in4.txt -meta ".data_folder().$name."_in4_meta.txt -type cytoscan_hd -family TEST -out_map $map -out_ped $ped -out_fam $fam");
check_file($fam, data_folder().$name."_out4_fam.txt");
check_file($ped, data_folder().$name."_out4_ped.txt");
check_file($map, data_folder().$name."_out4_map.txt");

end_test();

?>