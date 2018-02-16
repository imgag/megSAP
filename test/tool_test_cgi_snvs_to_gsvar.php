<?php

require_once("framework.php");

$name = "cgi_snvs_to_gsvar";

start_test($name);

//1. test: input file without known CGI cancer acronym translation
$gsvar_input1 = data_folder() . "cgi_snvs_to_gsvar_in1.GSvar";
$gsvar_ref1 = data_folder() . "cgi_snvs_to_gsvar_ref1.GSvar";
$cgi_mutation_file_input1 = data_folder() . "cgi_snvs_to_gsvar_in1_cgi_mutation_analysis.tsv";
$gsvar_output1 = output_folder() . "cgi_snvs_to_gsvar_out1.GSvar"; 

check_exec("php ".src_folder()."/NGS/cgi_snvs_to_gsvar.php -gsvar_in $gsvar_input1 -cgi_snv_in $cgi_mutation_file_input1 -out $gsvar_output1");
check_file($gsvar_output1,$gsvar_ref1,true);

//2. test: input file with known CGI cancer acronym translation
$gsvar_input2 = data_folder() . "cgi_snvs_to_gsvar_in2.GSvar";
$gsvar_ref2 = data_folder() . "cgi_snvs_to_gsvar_ref2.GSvar";
$cgi_mutation_file_input2 = data_folder() . "cgi_snvs_to_gsvar_in2_cgi_mutation_analysis.tsv";
$gsvar_output2 = output_folder() . "cgi_snvs_to_gsvar_out2.GSvar"; 
check_exec("php ".src_folder()."/NGS/cgi_snvs_to_gsvar.php -gsvar_in $gsvar_input2 -cgi_snv_in $cgi_mutation_file_input2 -out $gsvar_output2");
check_file($gsvar_output2,$gsvar_ref2,true);



?>
