<?php
require_once("framework.php");

$name = "converter_manta2tsv";
$script = realpath(src_folder())."/Tools/".$name.".php";
start_test($name);

$input1 = data_folder()."/converter_manta2tsv_in1.vcf.gz";
$ref1 = data_folder()."/converter_manta2tsv_ref1.tsv";
$out1 = output_folder()."/converter_manta2tsv_out1.tsv";

check_exec("php {$script} -in {$input1} -out {$out1} -tumor_id DX000001_01");
check_file($out1,$ref1);

end_test();
?>