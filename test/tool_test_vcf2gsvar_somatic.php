<?php

require_once("framework.php");

$name = "vcf2gsvar_somatic";
start_test($name);

// standard strelka
$out_file1 = output_folder().$name."_out1.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in1.vcf -out $out_file1 -t_col Tumor -n_col Normal --log ".output_folder().$name."_out1.log");
check_file($out_file1, data_folder().$name."_out1.GSvar");

// standard ion torrent
$out_file1 = output_folder().$name."_out2.GSvar";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in2.vcf.gz -out $out_file1 -t_col Tumor --log ".output_folder().$name."_out2.log");
check_file($out_file1, data_folder().$name."_out2.GSvar");

end_test();

?>
