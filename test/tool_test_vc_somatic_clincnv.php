<?php
require_once("framework.php");

$name = "vc_somatic_clincnv";
start_test($name);
//prepare input data_folder
$cov_folder = output_folder()."/cov_folder/";
$cohort_folder1 = output_folder()."/cohorts1/";
$cohort_folder2 = output_folder()."/cohorts2/";
exec2("mkdir -p $cov_folder");
exec2("rm -rf $cov_folder/*");
exec2("tar xzf ".data_folder().$name."_cov_files.tar.gz -C $cov_folder");
create_directory($cohort_folder1);
create_directory($cohort_folder2);

//test output folders for
$cov_folder_n = $cov_folder ."cov-normal/";
$cov_folder_t = $cov_folder ."cov-tumor/";
//file with tumor-normal pairs
$t_n_pair_file = $cov_folder_t."/list_tid-nid.csv";

//system files
$system_file = $cov_folder . "processing_system.ini";
$bed_file = $cov_folder ."target_region.bed";

//output file which will be compared to reference
$out_file1 = output_folder().$name."_out1.tsv";

$t_cov = "{$cov_folder_t}/DX000002_01.cov";
$n_cov = "{$cov_folder_n}/DX000002_02.cov";

//1st test: without off-target regions
check_exec("php ".src_folder()."/NGS/{$name}.php -t_id DX000002_01 -n_id DX000002_02 -cov_folder_n $cov_folder_n -cov_folder_t $cov_folder_t -cohort_folder $cohort_folder1 -cov_pairs $t_n_pair_file -out $out_file1 -t_cov $t_cov -n_cov $n_cov -bed $bed_file -system $system_file");
check_file($out_file1, data_folder().$name."_out1.tsv");

//2nd test: include off target regions
$bed_file_off = $cov_folder . "off_target.bed";
$out_file2 = output_folder().$name."_out2.tsv";
$cov_folder_n_off = $cov_folder ."cov-normal_off_target/";
$cov_folder_t_off = $cov_folder ."cov-tumor_off_target/";
$t_cov_off = $cov_folder_t_off . "DX000002_01.cov";
$n_cov_off = $cov_folder_n_off . "DX000002_02.cov";

check_exec("php ".src_folder()."/NGS/{$name}.php -t_id DX000002_01 -n_id DX000002_02 -cov_folder_n $cov_folder_n -cov_folder_t $cov_folder_t -cohort_folder $cohort_folder2 -cov_pairs $t_n_pair_file -out $out_file2 -t_cov $t_cov -n_cov $n_cov -bed $bed_file -system $system_file -t_cov_off $t_cov_off -n_cov_off $n_cov_off -cov_folder_n_off $cov_folder_n_off -cov_folder_t_off $cov_folder_t_off -bed_off $bed_file_off");

check_file($out_file2, data_folder().$name."_out2.tsv");

end_test();

?>