<?php
require_once("framework.php");

$name = "hla_genotyper";
start_test($name);

$bam = data_folder()."hla_genotyper_in.bam";

//test
$out_res = output_folder().$name."_out1.txt";
$out_dose = output_folder().$name."_out1.dose";
check_exec("php ".src_folder()."/NGS/{$name}.php -bam {$bam} -out {$out_res} -out_dose {$out_dose} -ethnicity EUR -name testSample_01 -build GRCh37 -genes A");
check_file($out_res, data_folder()."/{$name}_ref1.txt");
check_file($out_dose, data_folder()."/{$name}_ref1.dose");

end_test();

?>