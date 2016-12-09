<?php

require_once("framework.php");

$name = "plink_homocygosity";
start_test($name);

$out_file = output_folder().$name.".txt";
check_exec("php ".src_folder()."/Chips/".$name.".php -in_map ".data_folder().$name."_in_map.txt -in_ped ".data_folder().$name."_in_ped.txt -out $out_file -w_snp 60 -w_kb 2000 -w_het 8 -w_miss 10 -w_thres 0.5 -s_snp 10 -s_kb 1000 -s_dens 100 -s_gap 5000");
check_file($out_file, data_folder().$name."_out.txt");

end_test();

?>