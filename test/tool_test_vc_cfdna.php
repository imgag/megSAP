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
check_exec("php ".src_folder()."/NGS/{$name}.php -bam ".data_folder().$name."_in1.bam -target ".data_folder().$name."_in1_roi.bed -build GRCh38 -folder ".output_folder()." -monitoring_vcf ".data_folder().$name."_in1_monitoring.vcf");
//remove file date
exec2("sed -i '/^##fileDate/d' ".$vcf);
exec2("sed -i '/^##fileDate/d' ".$vcf_hq);
check_file($vcf, data_folder().$name."_out1.vcf");
check_file($vcf_hq, data_folder().$name."_out1_hq.vcf");
check_file($tsv_id, data_folder().$name."_out1_ID.tsv");
check_file($tsv_mon, data_folder().$name."_out1_monitoring.tsv");
end_test();

?>


