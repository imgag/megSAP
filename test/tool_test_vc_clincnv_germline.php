<?php
require_once("framework.php");

$name = "vc_clincnv_germline";

start_test($name);

//prepare input data_folder
$tmp_folder = output_folder()."/vc_clincnv_germline_data/";
exec2("rm -rf $tmp_folder/*");
exec2("tar xzf ".data_folder()."{$name}_data.tgz -C ".output_folder());
$bed = $tmp_folder."/target_region.bed";
$cov_folder = $tmp_folder."/coverage/";

//test
$out_file1 = output_folder().$name."_out1.tsv";
check_exec("php ".src_folder()."/NGS/{$name}.php -cov {$cov_folder}/DX000018_02.cov -cov_folder {$cov_folder} -cov_min 20 -bed {$bed} -out {$out_file1}");
check_file($out_file1, data_folder().$name."_out1.tsv");
check_file(substr($out_file1,0,-4).".seg", data_folder().$name."_out1.seg");
check_file(substr($out_file1,0,-4)."_cnvs.seg",data_folder().$name."_out2.seg");

end_test();

?>
