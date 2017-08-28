<?php

require_once("framework.php");

$name = "vc_manta";
start_test($name);

$t_bam = data_folder().$name."_tu_in.bam";
$n_bam = data_folder().$name."_no_in.bam";

//somatic mode
$out_file = output_folder().$name."_out_somatic.vcf.gz";
$out_file_indel = output_folder().$name."_out_somatic_smallindels.vcf.gz";

check_exec("php ".src_folder()."/NGS/vc_manta.php -t_bam $t_bam -bam $n_bam -out $out_file -smallIndels $out_file_indel -config_preset high_sensitivity -regions chr12:1-10000 --log ".output_folder().$name."_somatic.log");
check_file($out_file, data_folder().$name."_out_somatic.vcf.gz");
check_file($out_file_indel, data_folder().$name."_out_somatic_smallindels.vcf.gz");

//single BAM
$out_file_single = output_folder().$name."_out_single.vcf.gz";
$out_file_single_indel = output_folder().$name."_out_single_smallindels.vcf.gz";
check_exec("php ".src_folder()."/NGS/vc_manta.php -bam $n_bam -out $out_file_single -smallIndels $out_file_single_indel -regions chr12 --log ".output_folder().$name."_single.log");
check_file($out_file_single, data_folder().$name."_out_single.vcf.gz");
check_file($out_file_single_indel, data_folder().$name."_out_single_smallindels.vcf.gz");

end_test();