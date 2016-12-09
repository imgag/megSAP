<?php

require_once("framework.php");
require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../src/Common/genomics.php");

$name = "indel_realign";
start_test($name);

//test (tests two regions in one run to reduce runtime. See BED file for region coordinates.)
$prefix = output_folder().$name;
$out_file1 = output_folder().$name."_out1.bam";
check_exec("php ".src_folder()."/NGS/{$name}.php -in ".data_folder().$name."_in1.bam -prefix $prefix -genome ".get_path("data_folder")."/genomes/GATK/hg19.fa --log ".output_folder().$name."_out1.log");
// Change to canonical output filename
rename("${prefix}.realign.bam", $out_file1);
rename("${prefix}.realign.bam.bai", "${out_file1}.bai");
check_file($out_file1, data_folder().$name."_out1.bam");

end_test();

?>