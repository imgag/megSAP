<?php

///@todo remove once the evaluation of strelka2 is finished (MS)

require_once("framework.php");

$name = "vc_strelka";
start_test($name);

$t_bam = data_folder().$name."_tu_in.bam";
$n_bam = data_folder().$name."_no_in.bam";

$out_file = output_folder().$name."_out.vcf.gz";
check_exec("php ".src_folder()."/NGS/vc_strelka.php -t_bam $t_bam -n_bam $n_bam -out $out_file -k -amplicon --log ".output_folder().$name."_1.log");
check_file($out_file, data_folder().$name."_out.vcf.gz");

$out_file2 = output_folder().$name."_out2.vcf.gz";
check_exec("php ".src_folder()."/NGS/vc_strelka.php -t_bam $t_bam -n_bam $n_bam -out $out_file2 -amplicon --log ".output_folder().$name."_2.log");
check_file($out_file2, data_folder().$name."_out2.vcf.gz");

end_test();
?>
