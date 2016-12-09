<?php

require_once("framework.php");

$name = "chip_ngs_compare";
start_test($name);

//affymetrix_6.0
$out_file = output_folder().$name."_out.txt";
check_exec("php ".src_folder()."/Chips/".$name.".php -chip ".data_folder().$name."_in_chip.tsv -ngs ".data_folder().$name."_in_ngs.tsv -out $out_file -type affymetrix_6.0");
check_file($out_file, data_folder().$name."_out.txt");

end_test();

?>