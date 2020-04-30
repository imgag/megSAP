<?php

require_once("framework.php");

$name = "vc_cnvhunter";
start_test($name);

//prepare input data
$cov_nor = output_folder()."GS140794_02.cov";
copy(data_folder().$name."_in_nor.cov", $cov_nor);

$cov_folder = output_folder()."/cov_folder/";
exec2("mkdir -p $cov_folder");
exec2("rm -rf $cov_folder/*");
exec2("tar -m -xzf ".data_folder().$name."_cov_files.tgz -C $cov_folder");

//test germline (sample is supposed to be remission, but contains 1230 exon deletion chr6:28478466-33424368)
$out_file1 = output_folder().$name."_out1.tsv";
$out_file2seg = output_folder().$name."_out1.seg";
check_exec("php ".src_folder()."/NGS/{$name}.php -n 20  -cov $cov_nor --log ".output_folder().$name."_out2.log -system ".data_folder().$name."_system.ini -out $out_file1 -seg GS140794_02 -cov_folder $cov_folder");
// results depend on the availability of the NGSD 
if (db_is_enabled("NGSD"))
{
	check_file($out_file1, data_folder().$name."_out1.tsv");
}
else
{
	check_file($out_file1, data_folder().$name."_out1_noNGSD.tsv");
}
check_file($out_file2seg, data_folder().$name."_out1.seg");

end_test();

?>
