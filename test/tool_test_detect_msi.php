<?php

require_once("framework.php");

$name = "detect_msi";
start_test($name);

//delete status files
$in_t = data_folder()."detect_msi_tu_in.bam";
$in_n = data_folder()."detect_msi_no_in.bam";
$msi_ref = data_folder()."detect_msi_sites.list";
$out = output_folder().$name."_out1.tsv";
check_exec("php ".src_folder()."/Tools/{$name}.php -n_bam $in_n -t_bam $in_t -msi_ref $msi_ref -out $out");
check_file($out, data_folder().$name."_out1.tsv", true);
check(file_exists(output_folder().$name."_out1.tsv_somatic"), false);

//keep status files
$in_t = data_folder()."detect_msi_tu_in.bam";
$in_n = data_folder()."detect_msi_no_in.bam";
$msi_ref = data_folder()."detect_msi_sites.list";
$out = output_folder().$name."_out2.tsv";
check_exec("php ".src_folder()."/Tools/{$name}.php -n_bam $in_n -t_bam $in_t -msi_ref $msi_ref -out $out -keep_status_files");
check_file($out, data_folder().$name."_out2.tsv", true);
check_file($out."_dis", data_folder().$name."_out2.tsv_dis", true);
check_file($out."_all", data_folder().$name."_out2.tsv_all", true);
check_file($out."_unstable", data_folder().$name."_out2.tsv_unstable", true);

end_test();
?>