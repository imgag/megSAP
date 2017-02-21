<?php

require_once("framework.php");
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/genomics.php");

$name = "indel_realign";
start_test($name);

//test (tests two regions in one run to reduce runtime. See BED file for region coordinates.)
$out_file1 = output_folder().$name."_out1.bam";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in1.bam -out $out_file1 --log ".output_folder().$name."_out1.log");
check_file($out_file1, data_folder().$name."_out1.bam");

end_test();

?>