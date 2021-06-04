<?php

require_once("framework.php");

$name = "vc_cfdna";
start_test($name);

########################## clear output folder ##############
exec2("rm -f ".output_folder()."/*");

########################## VC ##################
$out_folder = output_folder();
$vcf = output_folder()."/${name}_in1.vcf";
$vcf_hq = output_folder()."/${name}_in1_hq.vcf";
$tsv_id = output_folder()."/${name}_in1_ID.tsv";
$tsv_mon = output_folder()."/${name}_in1_monitoring.tsv";
//print "php ".src_folder()."/NGS/{$name}.php -bam ".data_folder().$name."_in1.bam -target ".data_folder().$name."_in1_roi.bed -build GRCh37 -folder ".output_folder()." -monitoring_bed ".data_folder().$name."_in1_monitoring.bed --log ".output_folder().$name."_out1.log";
check_exec("php ".src_folder()."/NGS/{$name}.php -bam ".data_folder().$name."_in1.bam -target ".data_folder().$name."_in1_roi.bed -build GRCh37 -folder ".output_folder()." -monitoring_bed ".data_folder().$name."_in1_monitoring.bed");
check_file($vcf, data_folder().$name."_out1.vcf");
check_file($vcf_hq, data_folder().$name."_out1_hq.vcf");
check_file($tsv_id, data_folder().$name."_out1_ID.tsv");
check_file($tsv_mon, data_folder().$name."_out1_monitoring.tsv");
end_test();

?>


