<?php

require_once("framework.php");

$name = "an_somatic_cnvs";

start_test($name);

//NCG gene annotation
$cnv_input_1 = data_folder()."/an_somatic_cnvs_in1_cnvs.tsv"; 
$cnv_ref_1 = data_folder()."/an_somatic_cnvs_ref1_cnvs.tsv";
$cnv_output_1 = output_folder()."/an_somatic_cnvs_out1_cnvs.tsv";

check_exec("php ".src_folder()."/NGS/an_somatic_cnvs.php -cnv_in $cnv_input_1  -include_ncg -out $cnv_output_1 ");
check_file($cnv_output_1,$cnv_ref_1,true);

//annotate RNA data
$cnv_input_2 = data_folder()."/an_somatic_cnvs_in2_cnvs.tsv"; 
$rna_counts = data_folder()."/an_somatic_cnvs_rna_counts.tsv";
$cnv_output_2 = output_folder()."/an_somatic_cnvs_out2_cnvs.tsv";
$cnv_ref_2 = data_folder()."/an_somatic_cnvs_ref2_cnvs.tsv";
check_exec("php ".src_folder()."/NGS/an_somatic_cnvs.php -cnv_in $cnv_input_2 -rna_counts $rna_counts -out $cnv_output_2 -rna_id RX01_01 -rna_ref_tissue colon");
check_file($cnv_output_2,$cnv_ref_2,false);

//5th case: annotate cytobands
$cnv_input_3 = $cnv_input_2;
$cnv_output_3 = output_folder()."/an_somatic_cnvs_out3_cnvs.tsv";
$cnv_ref_3 = data_folder()."/an_somatic_cnvs_ref3_cnvs.tsv";
check_exec("php ".src_folder()."/NGS/an_somatic_cnvs.php -cnv_in $cnv_input_3 -out $cnv_output_3 -include_cytoband");
check_file($cnv_output_3,$cnv_ref_3,false);
end_test();
?>