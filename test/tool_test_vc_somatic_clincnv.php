<?php
require_once("framework.php");

$name = "vc_somatic_clincnv";
start_test($name);
//prepare input data_folder
$cov_folder = output_folder()."/cov_folder/";
$cohort_folder = output_folder()."/cohorts/";
exec2("mkdir -p $cov_folder");
exec2("rm -rf $cov_folder/*");
exec2("tar xzf ".data_folder().$name."_cov_files.tar.gz -C $cov_folder");
create_directory($cohort_folder);

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

check_exec("php ".src_folder()."/NGS/{$name}.php -t_id DX000002_01 -n_id DX000002_02 -cov_folder_n $cov_folder_n -cov_folder_t $cov_folder_t -cohort_folder $cohort_folder -cov_pairs $t_n_pair_file -out $out_file1 -t_cov $t_cov -n_cov $n_cov -bed $bed_file -system $system_file");
check_file($out_file1, data_folder().$name."_out1.tsv");

//check whether each tumor-normal pair from $t_n_pair_file was created
foreach(file($t_n_pair_file) as $line)
{
	if(strpos($line,"#") !== false) continue;
	list($t_id,$n_id) = explode(',',trim($line));
	$cnv_file = $cohort_folder . "/somatic/" .$t_id . "-" . $n_id . "/CNAs.txt";
	check_file_exists($cnv_file);
}
end_test();
?>