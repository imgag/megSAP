<?php

require_once("framework.php");

$name = "cgi_snvs_to_gsvar";

start_test($name);

$gsvar_input = data_folder() . "cgi_snvs_to_gsvar_in1.GSvar";
$gsvar_ref = data_folder() . "cgi_snvs_to_gsvar_ref1.GSvar";
$cgi_mutation_file_input = data_folder() . "cgi_snvs_to_gsvar_in1_cgi_mutation_analysis.tsv";
$gsvar_output = output_folder() . "cgi_snvs_to_gsvar_out1.GSvar"; 

check_exec("php ".src_folder()."/NGS/cgi_snvs_to_gsvar.php -gsvar_in $gsvar_input -cgi_snv_in $cgi_mutation_file_input -out $gsvar_output");
check_file($gsvar_output,$gsvar_ref,true);



?>
