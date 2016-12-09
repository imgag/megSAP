<?php

require_once("framework.php");
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/genomics.php");

$name = "indel_realign_abra";
start_test($name);

//test (tests two regions in one run to reduce runtime. See indel_realign_abra test BED file for region coordinates.)
$out_file1 = output_folder().$name."_out1.bam";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in1.bam -roi ".data_folder().$name."_roi.bed -out $out_file1 -threads 1 -mer 0.1 -build hg19 --log ".output_folder().$name."_out1.log");
check_file($out_file1, data_folder().$name."_out1.bam");

end_test();

?>