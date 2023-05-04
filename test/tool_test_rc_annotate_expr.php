<?php

require_once("framework.php");

$name = "rc_annotate_expr";
start_test($name);

//test1
$out = output_folder().$name."_out.tsv";
$out_cohort = output_folder().$name."_out_cohort.tsv";
$out_stats = output_folder().$name."_out_stats.tsv";
$out_corr = output_folder().$name."_out_corr.txt";
$in = data_folder().$name."_in1.tsv";
$in_files = data_folder().$name."_in_files.tsv";

$out_hpa = output_folder().$name."_out_hpa.tsv";
$in_hpa = data_folder().$name."_in_hpa.tsv";
$out_hpa_corr = output_folder().$name."_out_corr_hpa.txt";

//test 1: annotate with cohort
check_exec("php ".src_folder()."/NGS/".$name.".php -name Sample1 -in $in -in_files $in_files -out $out -cohort $out_cohort -stats $out_stats -corr $out_corr --log ".output_folder().$name."_1.log");
check_exec("numdiff -a 0.1 " . $out . " " . data_folder().$name."_out.tsv");
check_exec("numdiff -a 0.1 " . $out_cohort . " " . data_folder().$name."_out_cohort.tsv");
check_exec("numdiff -a 0.1 " . $out_stats . " " . data_folder().$name."_out_stats.tsv");
check_exec("numdiff -a 0.1 " . $out_corr . " " . data_folder().$name."_out_corr.txt");

//test 2: annotate with cohort and HPA reference
check_exec("php ".src_folder()."/NGS/".$name.".php -name Sample1 -in $in -in_files $in_files -out $out_hpa -cohort $out_cohort -stats $out_stats -corr $out_corr -hpa_corr $out_hpa_corr -hpa_ref $in_hpa -hpa_tissue colon");
check_exec("numdiff -a 0.1 " . $out_hpa . " " . data_folder().$name."_out_hpa.tsv");
check_exec("numdiff -a 0.1 " . $out_cohort . " " . data_folder().$name."_out_cohort.tsv");
check_exec("numdiff -a 0.1 " . $out_stats . " " . data_folder().$name."_out_stats.tsv");
check_exec("numdiff -a 0.1 " . $out_corr . " " . data_folder().$name."_out_corr.txt");
check_exec("numdiff -a 0.1 " . $out_hpa_corr . " " . data_folder().$name."_out_corr_hpa.txt");

end_test();

?>