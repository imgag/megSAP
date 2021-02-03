<?php
require_once("framework.php");

$name = "vc_clincnv_germline";
$ngsbits = get_path("ngs-bits");

start_test($name);

//test germline 1
//prepare input data_folder
$tmp_folder = output_folder()."/vc_clincnv_germline_data1/";
exec2("rm -rf $tmp_folder/*");
exec2("tar -m -xzf ".data_folder()."{$name}_data1.tar.gz -C ".output_folder());
$bed = $tmp_folder."/target_region.bed";
$cov_folder = $tmp_folder."/coverage/";

//test1
$out_file1 = output_folder().$name."test1_out.tsv";
$log_file1 = output_folder().$name."test1_out.log";
check_exec("php ".src_folder()."/NGS/{$name}.php -cov {$cov_folder}/DX000018_02.cov -cov_folder {$cov_folder} -cov_min 20 -bed {$bed} -out {$out_file1} --log {$log_file1}");
check_file($out_file1, data_folder().$name."_test1_out.tsv");
check_file(substr($out_file1,0,-4).".seg", data_folder().$name."_test1_out1.seg");
check_file(substr($out_file1,0,-4)."_cnvs.seg",data_folder().$name."_test1_out2.seg");

//test germline 2 (old CnvHunter test)
//prepare input data_folder
$cov_folder = output_folder()."/vc_clincnv_germline_data2/";
exec2("mkdir -p $cov_folder");
exec2("rm -rf $cov_folder/*");
exec2("tar -m -xzf ".data_folder()."{$name}_data2.tar.gz -C ".output_folder());
$bed = $cov_folder."/ssHAEv5_2016_07_11_annotated.bed";
//test2
$out_file2 = output_folder().$name."_test2_out.tsv";
$log_file2 = output_folder().$name."_test2_out.log";
check_exec("php ".src_folder()."/NGS/{$name}.php -cov {$cov_folder}/GS140794_02.cov -cov_folder {$cov_folder} -cov_min 20 -max_cnvs 200 -bed {$bed} -out {$out_file2} --log {$log_file2}");
check_file($out_file2, data_folder().$name."_test2_out.tsv");
check_file(substr($out_file2,0,-4).".seg", data_folder().$name."_test2_out1.seg");
check_file(substr($out_file2,0,-4)."_cnvs.seg",data_folder().$name."_test2_out2.seg");

//test clincnv for somatic single sample (vc_clincnv_germline is used in this case)
//prepare input data_folder
$tmp_folder = output_folder()."/vc_clincnv_tumor_data/";
exec2("mkdir -p $tmp_folder");
exec2("rm -rf $tmp_folder/*");
exec2("tar -xzf ".data_folder()."vc_clincnv_somatic_files.tar.gz -C {$tmp_folder}");
$bed_in = $tmp_folder."/target_region.bed";
$cov_folder = $tmp_folder."cov-tumor";
$bed = $tmp_folder."/target_region_annotated.bed";
$pipeline = [
        ["{$ngsbits}BedAnnotateGC", "-in {$bed_in} -ref /mnt/share/data/genomes/GRCh37.fa"],
        ["{$ngsbits}BedAnnotateGenes", "-out {$bed}"],
    ];
//off target
$bed_off = $tmp_folder . "off_target.bed";
$cov_off = $tmp_folder . "cov-tumor_off_target/DX000015_01.cov";
$parser = new ToolBase("tool_test_vc_clincnv_germline", "Pipeline for Bed Annotation.");
$parser->execPipeline($pipeline, "creating annotated BED file for ClinCNV");
$cov_folder = $tmp_folder."/cov-tumor";
//test
$out_file3 = output_folder().$name."_tumor_out.tsv";
$log_file3 = output_folder().$name."_tumor_out.log";
check_exec("php ".src_folder()."/NGS/{$name}.php -cov {$cov_folder}/DX000015_01.cov -cov_folder {$cov_folder} -cov_max 200 -max_cnvs 200 -bed {$bed} -bed_off {$bed_off} -cov_off {$cov_off} -out {$out_file3} --log {$log_file3} -tumor_only");
check_file($out_file3, data_folder().$name."_tumor_out.tsv");
check_file(substr($out_file3,0,-4).".seg", data_folder().$name."_tumor_out1.seg");
check_file(substr($out_file3,0,-4)."_cnvs.seg",data_folder().$name."_tumor_out2.seg");

end_test();

?>
