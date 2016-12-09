<?php

require_once("framework.php");

$name = "create_samplesheet";
start_test($name);

$out_file = output_folder().$name."_out1.csv";
check_exec("php ".src_folder()."/NGS/".$name.".php -in ".data_folder().$name."_in1.tsv -out $out_file -fcid C3NPWACXX -run 304");
check_file($out_file, data_folder().$name."_out1.csv");

end_test();

?>