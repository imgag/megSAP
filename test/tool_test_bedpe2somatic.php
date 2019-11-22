<?php

require_once("framework.php");

$name = "bedpe2somatic";
start_test($name);

$out = output_folder().$name."_out1.bedpe";
$in = data_folder().$name."_in1.bedpe";
$ref = data_folder().$name."_ref1.bedpe";
check_exec("php ".src_folder()."/Tools/bedpe2somatic.php -in $in -out $out -tid DX000002_01 -nid DX000001_01");
check_file($out, $ref, true);

$out = output_folder().$name."_out2.bedpe";
$in = data_folder().$name."_in2.bedpe";
$ref = data_folder().$name."_ref2.bedpe";
check_exec("php ".src_folder()."/Tools/bedpe2somatic.php -in $in -out $out -tid DX000002_01 -nid DX000001_01");
check_file($out, $ref, true);


end_test();
?>