<?php

require_once("framework.php");

$name = "vc_strelka2";
start_test($name);

$t_bam = data_folder().$name."_tu_in.bam";
$n_bam = data_folder().$name."_no_in.bam";

$out_file = output_folder().$name."_out.vcf.gz";
check_exec("php ".src_folder()."/NGS/vc_strelka2.php -t_bam $t_bam -n_bam $n_bam -out $out_file -k -smallIndels ".data_folder()."/vc_strelka2_smallIndels_in.vcf.gz --log ".output_folder().$name."_1.log");
check_file($out_file, data_folder().$name."_out.vcf.gz");

end_test();
?>
