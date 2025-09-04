<?php

require_once("framework.php");
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/genomics.php");

$name = "rc_normalize_cohort_exons";

start_test($name);

$in_cohort1 = data_folder()."{$name}_in_cohort1.tsv.gz";
$in_raw_counts1 = data_folder()."{$name}_in_sampleraw1.tsv";
$in_norm_genes1 = data_folder()."{$name}_in_norm_genes1.tsv";
$out_log1 = output_folder().$name."_log.txt";
$out_cohort1 = output_folder().$name."_out_cohort1.tsv.gz";
$out_sample1 = output_folder().$name."_out_sample1.tsv";
check_exec("php ".src_folder()."/Tools/{$name}.php -ps_name Sample1 -in_cohort {$in_cohort1} -normalized_genes {$in_norm_genes1} -uncorrected_counts {$in_raw_counts1} -out_cohort {$out_cohort1} -out_sample {$out_sample1} --log {$out_log1}");
check_file($out_cohort1, data_folder().$name."_out_cohort1.tsv.gz");
check_file($out_sample1, data_folder().$name."_out_sample1.tsv");

end_test();

?>