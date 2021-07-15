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
check_exec("php ".src_folder()."/NGS/".$name.".php -name Sample1 -in $in -in_files $in_files -out $out -cohort $out_cohort -stats $out_stats -corr $out_corr");
check_file($out, data_folder().$name."_out.tsv", true, "1e-6");
check_file($out_cohort, data_folder().$name."_out_cohort.tsv", true, "1e-6");
check_file($out_stats, data_folder().$name."_out_stats.tsv", true, "1e-6");
check_file($out_corr, data_folder().$name."_out_corr.txt", true, "1e-6");

//test 2: annotate with cohort and HPA reference
check_exec("php ".src_folder()."/NGS/".$name.".php -name Sample1 -in $in -in_files $in_files -out $out_hpa -cohort $out_cohort -stats $out_stats -corr $out_corr -hpa_corr $out_hpa_corr -hpa_ref $in_hpa -hpa_tissue colon");
check_file($out, data_folder().$name."_out_hpa.tsv", true, "1e-6");
check_file($out_cohort, data_folder().$name."_out_cohort.tsv", true, "1e-6");
check_file($out_stats, data_folder().$name."_out_stats.tsv", true, "1e-6");
check_file($out_corr, data_folder().$name."_out_corr.txt", true, "1e-6");
check_file($out_hpa_corr, data_folder().$name."_out_corr_hpa.txt", true, "1e-6");

end_test();

?>