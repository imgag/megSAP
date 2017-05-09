<?php

require_once("framework.php");
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/genomics.php");

$name = "indel_realign_abra";
start_test($name);

//test with ROI (tests two regions in one run to reduce runtime. See indel_realign_abra test BED file for region coordinates.)
$out_file1 = output_folder().$name."_out1.bam";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in1.bam -roi ".data_folder().$name."_roi.bed -out $out_file1 --log ".output_folder().$name."_out1.log");
check_file($out_file1, data_folder().$name."_out1.bam");

//test without ROI
$out_file2 = output_folder().$name."_out2.bam";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in1.bam -out $out_file2 --log ".output_folder().$name."_out2.log");
check_file($out_file2, data_folder().$name."_out2.bam");

end_test();

?>